/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/PDG.h"
#include "CepGenProcesses/DiffVM.h"

namespace CepGen {
  namespace Process {
    DiffVM::DiffVM()
        : GenericProcess("diffvm", "Diffractive vector meson production"),
          ifragp_(BeamMode::Elastic),
          ifragv_(BeamMode::Elastic),
          igammd_(PhotonMode::Fixed),
          bmin_(0.),
          dmxv_(mp_),
          dmxp_(mp_) {}

    void DiffVM::beforeComputeWeight() {
      if (ifragp_ == BeamMode::Elastic) {
        const auto& w_limits = cuts_.cuts.central.mass_single;
        if (!w_limits.hasMin())
          throw CG_FATAL("DiffVM") << "You must specify a lower limit to W(gamma,p)!\n\t"
                                   << "Current limits: " << w_limits << ".";
        bmin_ = slp_.b0 + 4. * pom_.alpha1 * log(w_limits.min() / slp_.wb0);
      } else {
        if (ifragv_ == BeamMode::Elastic)
          bmin_ = slp_.b0 + 4. * pom_.alpha1 * log(slp_.amxb0 / slp_.wb0);
        else
          bmin_ = slp_.b0 + 4. * pom_.alpha1 * log(4. * pow(slp_.amxb0, 2) / (slp_.wb0 * sqs_));
      }
      bmin_ = std::max(bmin_, 0.5);
      CG_DEBUG("DiffVM") << "Minimum b slope: " << bmin_ << ".";
    }

    void DiffVM::addEventContent() {
      GenericProcess::setEventContent({{Particle::IncomingBeam1, {PDG::proton}},
                                       {Particle::IncomingBeam2, {PDG::proton}},
                                       {Particle::Parton1, {PDG::photon}},
                                       {Particle::Parton2, {PDG::pomeron}},
                                       {Particle::OutgoingBeam1, {PDG::proton}},
                                       {Particle::OutgoingBeam2, {PDG::proton}},
                                       {Particle::CentralSystem, {PDG::Upsilon1S}}});
      event_->dump();
    }

    double DiffVM::computeWeight() {
      const double pout = 0., q2 = 0.;

      generatePhoton(x(0));

      switch (ifragp_) {
        case BeamMode::Elastic:
          break;
        default: {
          double dum = 0.;
          dmxp_ = outgoingPrimaryParticleMass(x(2), dum, true);
        } break;
      }
      switch (ifragv_) {
        case BeamMode::Elastic:
          break;
        default: {
          double dum = 0.;
          dmxv_ = outgoingPrimaryParticleMass(x(3), dum, false);
        } break;
      }

      //--- determine actual CM energy
      const Particle::Momentum p_pgam =
          event_->getOneByRole(Particle::IncomingBeam2).momentum() + event_->getOneByRole(Particle::Parton1).momentum();
      w2_ = p_pgam.energy2();
      const double w = sqrt(w2_);

      //--- return if generated masses are bigger than CM energy
      if (dmxp_ + dmxv_ > w - 0.1)
        return 0.;

      //--- calculate slope parameter b
      // generate t with e**(b*t) distribution

      double b = slp_.b0 + 4. * pom_.alpha1 * log(w / slp_.wb0);
      if (ifragp_ != BeamMode::Elastic)
        b -= 4. * pom_.alpha1m * log(dmxp_ / slp_.amxb0);
      if (ifragv_ != BeamMode::Elastic)
        b -= 4. * pom_.alpha1 * log(dmxv_ / slp_.amxb0);
      b = std::max(b, 0.5);

      const double t = computeT(x(1), b);
      CG_DEBUG_LOOP("DiffVM:weight") << "computed t=" << t << " GeV² for b=" << b << ".";

      //--- calculate actual minimal and maximal t for the generated masses
      // note that t here is positive!

      // Formula (E.5) from Review of Particle Properties 1992, p. III.50
      // 1: gamma, 2: p, 3: VM(+X), 4: p remnant
      // The formula for Pcm1 is altered to take the imaginary photon mass into account.

      const double inv_w = 1. / w;
      const double pcm1 = 0.5 * sqrt(pow(w2_ + q2 - mp2_, 2) + 4. * q2 * mp2_) * inv_w;
      const double pcm3 = 0.5 * sqrt((w2_ - pow(dmxv_ + dmxp_, 2)) * (w2_ - pow(dmxv_ - dmxp_, 2))) * inv_w;
      const double t_mean = 0.5 * ((-q2 - mp2_) * (dmxv_ * dmxv_ - dmxp_ * dmxp_) * inv_w * inv_w + w2_ + q2 - mp2_ -
                                   dmxv_ * dmxv_ - dmxp_ * dmxp_);
      const double t_min = t_mean - 2. * pcm1 * pcm3, t_max = t_mean + 2. * pcm1 * pcm3;
      if (t < t_min || t > t_max)
        return 0.;

      const double yhat = 0.25 * (t - t_min) / (pcm1 * pcm3);
      //      std::cout << bmin_ << std::endl;
      return bmin_ / b;

      //===== GENDIF

      const double ctheta = 1. - 2. * yhat, stheta = 2. * sqrt(yhat - yhat * yhat);

      //--- calculate 5-vectors of diffractive states in the CMS

      const double p_gamf = pout * ctheta / p_vm_.p();
      const double phi = 2. * M_PI * x(1);
      const Particle::Momentum pt(
          -cos(phi) * p_vm_.pz(), sin(phi) * p_vm_.pz(), -sin(phi) * p_vm_.py() + cos(phi) * p_vm_.px());
      const double ptf = pout * stheta / std::hypot(p_vm_.pz(), pt.pz());
      Particle::Momentum p_vm_cm = p_gamf * p_vm_ + ptf * pt;
      p_vm_cm.setMass(dmxv_);

      if (fabs(pout * pout - p_vm_cm.p2()) > 0.01 * pout * pout)
        CG_WARNING("DiffVM:weight") << "POUT <> |PCMVMX|\n\t"
                                    << "POUT: " << pout << ", PCMVMX = " << p_vm_cm << ".";

      Particle::Momentum p_px_cm = Particle::Momentum() - p_vm_cm;
      p_px_cm.setMass(dmxp_);

      //--- calculate momentum carried by the pomeron
      // the pomeron is thougt to be a quasireal particle emitted by the proton
      // and absorbed by the virtual vector meson

      //.....

      return 0.;
    }

    unsigned int DiffVM::numDimensions(const Kinematics::Mode&) const { return 4; }

    void DiffVM::fillKinematics() {}

    void DiffVM::generatePhoton(double x) {
      const Particle::Momentum p_ib = event_->getOneByRole(Particle::IncomingBeam1).momentum();
      //const Particle::Momentum p_ib2 = event_->getOneByRole( Particle::IncomingBeam2 ).momentum();
      Particle::Momentum p_pho, p_ob;
      if (igammd_ == PhotonMode::Fixed) {  // fixphot
        const double egamma = 3.;          // in steering card
        const double y = egamma / p_ib.energy();
        p_pho = y * p_ib;
        p_pho.setMass(-sqrt(fabs(p_ib.mass2() * y * y / (1. - y))));
        p_ob = p_ib - p_pho;
      }
      //else if ( igammd_ ==
      event_->getOneByRole(Particle::Parton1).setMomentum(p_pho);
      event_->getOneByRole(Particle::OutgoingBeam1).setMomentum(p_ob);
    }

    double DiffVM::outgoingPrimaryParticleMass(double x, double& y, bool treat) const {
      const auto& m_range = cuts_.cuts.remnants.mass_single;
      double m = 0.;
      if (fabs(pom_.epsilm) < 1.e-3) {
        //--- basic spectrum: 1/M^2
        const double lmin = 2. * log(m_range.min());
        const double delta = 2. * log(m_range.max() / m_range.min());
        m = sqrt(exp(x * delta + lmin));
      } else {
        //--- basic spectrum: 1/M^2(1+epsilon)
        const double m2min = pow(m_range.min(), -2. * pom_.epsilm);
        const double fact = pow(m_range.max(), -2. * pom_.epsilm) - m2min;
        m = sqrt(pow(fact * x + m2min, -1. / pom_.epsilm));
      }
      if (m < m_range.min()) {
        CG_ERROR("DiffVM:mass") << "M=" << m << " < MMIN=" << m_range.min() << ".";
        return m_range.min();
      }
      if (m > m_range.max()) {
        CG_ERROR("DiffVM:mass") << "M=" << m << " > MMAX=" << m_range.max() << ".";
        return m_range.max();
      }
      if (!treat)
        return m;

      const double m2 = m * m;
      //--- old version with enhancements in lower mass region
      if (m2 >= 4.)
        y = 1.;
      else if (m2 >= 3.1)
        y = 1.64 - 0.16 * m2;
      else if (m2 >= 2.65)
        y = m2 * (0.47 - 0.42 * pow(m2 - 2.65, 2));
      else if (m2 >= 2.25)
        y = m2 * (0.47 + 0.46 * pow(m2 - 2.65, 2));
      else if (m2 >= 2.02)
        y = m2 * (0.76 - 2.69 * pow(m2 - 2.02, 2));
      else if (m2 >= 1.72)
        y = m2 * (0.76 - 1.98 * pow(m2 - 2.02, 2));
      else
        y = 1.05 * (m2 - 1.165);

      if (1.60 * rand() < y)
        return m;

      /*//--- new version: 1/M_x^2 spectrum
      if ( m2 >= 2. ) y = 1.;
      else            y = 1.-1.1815*pow( m2-2., 2 );
      if ( rand() < y || ++i >= 100 )
        break;*/

      return 0.;
    }

    double DiffVM::computeT(double x, double b) const {
      const auto& t_range = cuts_.cuts.initial.q2;
      const double t_min = t_range.min(), t_max = t_range.max();

      double bloc = b;

      //--- generate spectrum by method of R. Lausen

      CG_DEBUG_LOOP("DiffVM:t") << "t range: " << t_range << ", b: " << b << ".";

      if (b < 0.1) {
        CG_ERROR("DiffVM:t") << "b = " << b << " < 0.1.";
        bloc = 0.1;
      }
      if (t_min >= t_max) {
        CG_ERROR("DiffVM:t") << "t range: " << t_range << "=> return " << t_min << ".";
        return t_min;
      }

      if (slp_.anexp < 1.) {
        // power law exponent is 0 or illegal => generate pure exp(bt) spectrum
        if (bloc * (t_max - t_min) >= 25.)  // method 1
          return t_min - log(x) / bloc;
        // method 2
        return t_min - log(1. - x * (1. - exp(bloc * (t_min - t_max)))) / bloc;
      }
      //--- new 16.5.07 BL:
      // Generate mixed exp(bt)/power law spectrum
      // d sigma/d t = exp (-n*ln(-bt/n+1)) = (-bt/n+1)^-n
      // Limit for small bt: exp (bt + c t^2) with c=b^2/2n
      // Limit for large bt>>n: t^-n
      const double c1 = pow(slp_.anexp + bloc * t_min, 1. - slp_.anexp);
      const double c0 = pow(slp_.anexp + bloc * t_max, 1. - slp_.anexp);
      const double t = -(slp_.anexp - pow(x * (c1 - c0) + c0, 1. / (1. - slp_.anexp))) / bloc;
      CG_DEBUG_LOOP("DiffVM:t") << "X=" << x << ", C0=" << c0 << ", C1=" << c1 << ", ANEXP=" << slp_.anexp
                                << ", BLOC=" << bloc << ", T=" << t;
      return t;
    }
  }  // namespace Process
}  // namespace CepGen
