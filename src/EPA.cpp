#include <iostream>

#include <CepGen/Core/Exception.h>
#include <CepGen/Physics/Constants.h>
#include "CepGenDiffVM/EPA.h"

namespace CepGen {
  const double EPA::ALPHARED = Constants::alphaEM * 0.5 * M_1_PI;

  EPA::EPA(const EPA::Mode& mode)
      : mode_(mode),
        y_range_({0., 1.}),
        dy_range_({0., 1.}),
        s_(0.),
        w12_(0.),
        elpr_(0.),
        eel_(0.),
        num_errors_({{0, 0}}),
        max_errors_({10, 500}),
        epa_max_(-1.) {}

  void EPA::init(const Momentum& pel, const Momentum& ppr, const Limits& q2_range, const Limits& w_range) {
    //--- define the kinematics ranges
    pel_ = pel;
    ppr_ = ppr;
    q2_range_ = q2_range;
    w_range_ = w_range;

    CG_INFO("EPA:init") << "EPA module initialised.\n\t"
                        << "mode: " << mode_ << "\n\t"
                        << "beams momenta:\n\t  " << pel_ << "\n\t  " << ppr_ << "\n\t"
                        << "W range: " << w_range_ << " GeV, "
                        << "Q² range: " << q2_range_ << " GeV².";

    //--- calculate CMS s=(P+K)²
    const double sqrt_s = CMEnergy(pel_, ppr_);
    s_ = sqrt_s * sqrt_s;
    elpr_ = 0.5 * (s_ - pel_.mass2() - ppr_.mass2());
    if (mode_ == Mode::transverse_longitudinal_pframe)
      //--- evaluate photon flux in proton rest frame:
      // set EEL to approx. 50TeV
      eel_ = elpr_ / ppr_.mass();
    else
      eel_ = pel_.energy();
    CG_DEBUG("EPA:init") << "S = " << s_ << ", EEL = " << eel_ << ".";
    const double w_min = w_range_.min(), w_max = w_range_.max();
    const double w_min2 = w_min * w_min;
    const double q2_min = q2_range_.min(), q2_max = q2_range_.max();
    w12_ = w_min2 - ppr_.mass2();

    //--- calculate Y - bounds from
    // ALI, A. et al. (1987): Heavy quark physics at HERA.
    // - Proc. HERA workshop, Hamburg 1987 (ed. R.D. PECCEI), 395-494.
    const double y_sqr = sqrt(pow(s_ - w12_, 2) - 4. * w12_ * pel_.mass2());
    dy_range_.max() = 0.5 * (s_ + w12_ + y_sqr) / (s_ + pel_.mass2());
    //---  use trick for quadratic equations; see
    // W.H. PRESS et al. (1988): Numerical Recipes in C. Cambridge
    // (Cambridge Univ. Press), p. 156.
    dy_range_.min() = std::max(y_range_.min(), w12_ / (dy_range_.max() * (s_ + pel_.mass2())));
    //--- calculate absolute maximum of y, irrespective of final state
    dy_range_.max() = std::min(
        y_range_.max(), std::min(s_ / (s_ + pel_.mass2()), 0.5 * (w_max * w_max - ppr_.mass2() + q2_max) / elpr_));
    CG_INFO("EPA:init") << "S / (S + DME**2) = " << (s_ / (s_ + pel_.mass2())) << "\n\t"
                        << "y range: " << y_range_ << "\n\t"
                        << "(W²(max)-mp²+Q²(max))/(2*elpr) = "
                        << (0.5 * (w_max * w_max - ppr_.mass2() + q2_max) / elpr_) << "\n\t"
                        << "dy range: " << dy_range_ << ".";
    const double dy_min = dy_range_.min(), dy_max = dy_range_.max();

    //--- set max. photon weight for efficient rejection plane
    const double qg2_min = std::max(q2_min, pel_.mass2() * pow(dy_min, 2) / (1. - dy_min));
    const double qg2_max = std::min(q2_max, dy_max * s_);
    Limits qg2_range(qg2_min, qg2_max);
    CG_DEBUG("EPA:init") << "QG2 range: " << qg2_range << ".";

    if (mode_ == Mode::wwa)  //--- WWA - approximation
      epa_max_ = ALPHARED * (4. * (1. - dy_min) + dy_min * dy_min);
    else {
      //--- full transversal spectrum (2) or full longitudinal and
      // transversal (3) spectrum
      const double eqe = qg2_min / eel_ / eel_;
      const double emqe2 = pow(dy_min - 0.25 * eqe, 2);
      const double emsqr =
          (pow(dy_min * elpr_, 2) + qg2_min * ppr_.mass2()) / (elpr_ * elpr_ + pel_.mass2() * ppr_.mass2());

      if (emsqr < 0.)
        CG_FATAL("EPA:init") << "problem with sqrt(emsqr): " << emsqr << " at EPAMAX determination.";

      epa_max_ = ALPHARED * dy_min * sqrt(emsqr) / (emqe2 + eqe);
      if (mode_ == Mode::transverse)
        epa_max_ *= (2. * (1. - dy_min) + emqe2 + eqe);
      else
        epa_max_ *= (4. * (1. - dy_min) + emqe2 + eqe);
    }
    epa_max_ *= log(dy_max / dy_min) * log(qg2_max / qg2_min);

    CG_DEBUG("EPA:init") << "maximal EPA: " << epa_max_;
  }

  EPA::Result EPA::operator()(double x_y, double x_q2) const {
    Result out;

    const double dy_min = dy_range_.min(), dy_max = dy_range_.max();

    //--- produce Y spect. ( 1/y weighted shape )
    const double y = dy_min * pow(dy_max / dy_min, x_y);
    //--- calculate actual Q2_min, Q2_max from Y
    const double gq2_min = std::max(pel_.mass2() * y * y / (1. - y), q2_range_.min());
    const double gq2_max = std::min(y * s_, q2_range_.max());
    const Limits gq2_range(gq2_min, gq2_max);
    //--- take Q2_cut from steering, if it is kinematicly reachable.
    if (!gq2_range.valid())
      return out;

    CG_DEBUG("EPA:get") << "Y, Q2min, Q2max: " << y << ", " << gq2_range << ".";

    double epa = 0., epa_t = 0., epa_l = 0., l_frac = 0.;

    //--- produce Q2 spect. (1/x weighted shape )
    const double q2 = gq2_min * pow(gq2_max / gq2_min, x_q2);

    //------------------------------------------------------------------
    // EPA - WWA spectrum
    //------------------------------------------------------------------

    //--- compute photon weight
    if (mode_ == Mode::wwa) {
      //--- WWA - approximation
      const double r = ALPHARED / (y * q2);
      epa_t = r * (2. * (1. - y) * (1. - pel_.mass2() * y * y / ((1. - y) * q2)) + y * y);
      epa_l = r * (2. * (1. - y));
      epa = epa_t + epa_l;
      l_frac = epa_l / epa;
    } else {
      //----------------------------------------------------------------
      // full transversal spectrum (2) or full longitudinal and
      // transversal (3) spectrum from:
      // * ABT, I. & J.R. SMITH (1992): MC Upgrades to Study
      //    Untagged Events. - H1-10/92-249.
      // See also:
      // * SMITH, J.R. (1992): An Experimentalist's Guide to Photon
      //    Flux Calculations. - H1-12/92-259.
      // * SMITH, J.R. (1993): Polarization Decomposition of Fluxes
      //    and Kinematics in ep Reactions. - H1-04/93-282.
      //----------------------------------------------------------------
      const double eqe = q2 / eel_ / eel_;
      const double emqe2 = pow(y - 0.25 * eqe, 2);
      const double emsqr = (pow(y * elpr_, 2) + q2 * ppr_.mass2()) / (elpr_ * elpr_ + pel_.mass2() * ppr_.mass2());

      if (emsqr < 0.) {
        CG_WARNING("EPA:get") << "problem with sqrt(emsqr) = " << emsqr << ": "
                              << "y, Q2 pair rejected";
        if (++num_errors_[0] > max_errors_[0])
          throw CG_FATAL("EPA:get") << "Errors threshold (" << max_errors_[0] << ") reached for "
                                    << "sqrt(emsqr) definition!\n\tTry WWA...";
      }

      const double r = ALPHARED / q2 * sqrt(emsqr) / (emqe2 + eqe);
      epa_t = r * (2. * (1. - y) + emqe2 + eqe);
      epa_l = (mode_ == Mode::transverse) ? 0. : r * (2. * (1. - y));  // longitudinal & transversal spectrum
    }

    epa = epa_t + epa_l;
    l_frac = epa_l / epa;

    //--- unweight MC
    double r = y * q2 * log(dy_max / dy_min) * log(gq2_max / gq2_min);
    const double w = sqrt(y * 2. * elpr_ - q2 + ppr_.mass2());
    //--- check if W_min < W < W_max, else reject photon
    if (!w_range_.passes(w))
      r = 0.;
    epa *= r;
    epa_t *= r;
    epa_l *= r;

    //--- update upper EPA bound
    if (epa > epa_max_) {
      if (epa > 1.1 * epa_max_)
        CG_WARNING("EPA:get") << "S: EPA > 1.1*max(EPA)!";
      else if (epa > 1.01 * epa_max_)
        CG_WARNING("EPA:get") << "W: EPA > 1.01*max(EPA)!";
      else
        CG_WARNING("EPA:get") << "I: EPA > max(EPA)!";
      epa_max_ = epa;
      CG_INFO("EPA:get") << "update of maximal weight: " << epa_max_ << ".";
    }

    CG_DEBUG_LOOP("EPA:get") << "Y: " << y << ", Q²: " << q2 << "\n\t"
                             << "EPA(T): " << epa_t << ", EPA(L): " << epa_l << "\n\t"
                             << "EPA: " << epa << ", long.fraction: " << l_frac << "\n\t"
                             << "GQ2 range: " << gq2_range << ".";

    //----> end rnd loop: rejection method
    if (rand() * 1. / RAND_MAX * epa_max_ > epa)
      return out;
    //--- continue with unweighted kin.

    //--- scattering angle of electron in LAB.:
    //  E.Lohrmann DESY HERA 83/08
    // x = Q2 / (y s)
    // E_sc = E_e(1-y) + E_p x y
    // cos t = [E_e(1-y) - E_p x y] / E_sc
    const double emy = pel_.energy() * (1. - y), exy = ppr_.energy() * q2 / s_;
    const double eesc = emy + exy;
    //const double cthe = ( emy-exy )/eesc, sthe = 2.*sqrt( emy-exy )/eesc;
    const double theta = acos((emy - exy) / eesc);

    //--- control scattering angle
    /*if ( fabs( cthe ) > 1. || fabs( sthe ) > 1. ) {
      CG_DEBUG_LOOP( "EPA:get" ) << "theta of electron: "
        << theta << ", "
        << "reject event for Y: " << y << ", Q2: " << q2 << ".";
      if ( ++num_errors_[1] > max_errors_[1] )
        throw CG_FATAL( "EPA:get" )
          << "Errors threshold (" << max_errors_[1] << ") reached for "
          << "theta definition!\n\t"
          << "theta: " << theta << ".";
      //--- new kinematics
      return out;
    }*/

    out.valid = true;
    const double phi = 2. * M_PI * rand() / RAND_MAX;

    //--- 5-vector of electron in LAB system
    const double pesc = -sqrt(eesc * eesc - pel_.mass2());
    out.ppe = Momentum::fromPThetaPhiE(pesc, theta, phi, eesc);

    //--- 5-vector of photon k - k'
    out.pph = pel_ - out.ppe;
    out.pph.setMass(sqrt(q2));

    //--- determine helicity of the photon
    out.heli = helicity(l_frac);
    return out;
  }

  short EPA::helicity(double longfr) const {
    if (rand() / RAND_MAX < longfr)
      return 0;
    if (rand() / 0.5)
      return +1;
    return -1;
  }

  std::ostream& operator<<(std::ostream& os, const EPA::Mode& mode) {
    switch (mode) {
      case EPA::Mode::wwa:
        return os << "WWA";
      case EPA::Mode::transverse:
        return os << "transverse only";
      case EPA::Mode::transverse_longitudinal:
        return os << "transverse+longitud.";
      case EPA::Mode::transverse_longitudinal_pframe:
        return os << "transverse+longitud.[p frame]";
    }
    return os;
  }
}  // namespace CepGen
