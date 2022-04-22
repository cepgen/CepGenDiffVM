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

#include <cassert>
#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/BreitWigner.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGenProcesses/DiffVM.h"

namespace cepgen {
  namespace proc {
    DiffVM::DiffVM(const ParametersList& params)
        : Process(params),
          vm_pdgid_(steer<ParticleProperties>("vmFlavour").pdgid),
          ifragp_(steerAs<int, BeamMode>("protonMode")),
          ifragv_(steerAs<int, BeamMode>("vmMode")),
          igammd_(steerAs<int, PhotonMode>("photonMode")),
          slp_(steer<ParametersList>("slopeParameters")),
          pom_(steer<ParametersList>("pomeronParameters")),
          vm_(steer<ParametersList>("vmParameters")),
          epa_calc_(steer<ParametersList>("epaParameters")) {}

    void DiffVM::prepareKinematics() {
      const Particle& ip2 = event_->oneWithRole(Particle::IncomingBeam2);
      mY2_ = ip2.mass2();

      const auto& w_limits = kin_.cuts().central.mass_single;
      const auto& q2_limits = kin_.cuts().initial.q2;

      //--- variables mapping
      defineVariable(
          phi_var_, Mapping::linear, {0., 2. * M_PI}, {0., 2. * M_PI}, "phi");  //FIXME extra factor 2*pi, to be checked
      defineVariable(t_var_, Mapping::linear, {0., 1.}, {0., 1.}, "Tvar");
      defineVariable(pho_var_, Mapping::linear, {0., 1.}, {0., 1.}, "PHOvar");
      defineVariable(difp_var_, Mapping::linear, {0., 1.}, {0., 1.}, "DPvar");

      if (igammd_ >= PhotonMode::WWA) {
        const Particle& ip1 = event_->oneWithRole(Particle::IncomingBeam1);
        epa_calc_.init(ip1.momentum(), ip2.momentum(), q2_limits, w_limits);
        defineVariable(wwa_var_, Mapping::linear, {0., 1.}, {0., 1.}, "WWAvar");
      }

      if (ifragp_ == BeamMode::Elastic) {
        if (!w_limits.hasMin())
          throw CG_FATAL("DiffVM") << "You must specify a lower limit to W(gamma,p)!\n\t"
                                   << "Current limits: " << w_limits << ".";
        bmin_ = slp_.b0 + 4. * pom_.alpha1 * log(w_limits.min() / slp_.wb0);
      } else {
        if (ifragv_ == BeamMode::Elastic)
          bmin_ = slp_.b0 + 4. * pom_.alpha1 * log(slp_.amxb0 / slp_.wb0);
        else {
          bmin_ = slp_.b0 + 4. * pom_.alpha1 * log(4. * pow(slp_.amxb0, 2) / (slp_.wb0 * sqs_));
          defineVariable(vm_var_, Mapping::linear, {0., 1.}, {0., 1.}, "VMvar");
        }
      }
      bmin_ = std::max(bmin_, 0.5);
      CG_DEBUG("DiffVM") << "Minimum b slope: " << bmin_ << ".";

      min_pho_energy_ = 0.25 * pow(w_limits.min(), 2) / ip2.momentum().p();
      max_s_ = pow(w_limits.max(), 2);

      if (vm_.lambda <= 0.)
        vm_.lambda = (*event_)[Particle::CentralSystem][0].mass();

      const double q2_min = q2_limits.min();
      prop_mx_ = std::max(1.,
                          vm_.xi * q2_min / (pow(vm_.lambda, 2) + vm_.xi * vm_.chi * q2_min) /
                              pow(1. + q2_min / pow(vm_.lambda, 2), vm_.eprop));

      const Particle& vm = (*event_)[Particle::CentralSystem][0];
      vm_mass_ = vm.mass(), vm_width_ = PDG::get().width(vm.pdgId());

      //--- mass range for VM generation
      double min_vm_mass = -1., max_vm_mass = -1.;
      const auto& invm_limits = kin_.cuts().central.mass_sum;
      if (invm_limits.valid()) {
        min_vm_mass = invm_limits.min();
        max_vm_mass = invm_limits.max();
      } else {
        min_vm_mass = vm_mass_ - 3. * vm_width_;
        max_vm_mass = vm_mass_ + 10. * vm_width_;
        if (vm.pdgId() == PDG::rho1450_0 || vm.pdgId() == PDG::rho1700_0)
          min_vm_mass = std::max(min_vm_mass, 1.2);
        else if (vm.pdgId() == PDG::h1380_1)
          min_vm_mass = std::max(min_vm_mass, 1.4);
      }
      vm_bw_.reset(new BreitWigner(vm_mass_, vm_width_, min_vm_mass, max_vm_mass));
    }

    void DiffVM::addEventContent() {
      GenericProcess::setEventContent({{Particle::IncomingBeam1, {PDG::electron}},
                                       {Particle::IncomingBeam2, {PDG::proton}},
                                       {Particle::Parton1, {PDG::photon}},
                                       {Particle::Parton2, {PDG::pomeron}},
                                       {Particle::OutgoingBeam1, {PDG::electron}},
                                       {Particle::OutgoingBeam2, {PDG::proton}},
                                       {Particle::CentralSystem, {vm_pdgid_}}});
    }

    double DiffVM::computeWeight() {
      //================================================================
      // GENGAM
      //================================================================

      if (!generatePhoton())
        return 0.;

      const double q2 = p_gam_.mass2();
      if (!kin_.cuts().initial.q2.contains(q2))
        return 0.;

      //--- determine gamma*-p energy
      p_cm_ = p_gam_ + event_->oneWithRole(Particle::IncomingBeam2).momentum();
      mB2_ = p_cm_.energy2();
      const double w = sqrt(mB2_);

      double weight = 1.;

      //--- determine weight of the virtual VM
      weight /= pow(1. + q2 / pow(vm_.lambda, 2), vm_.eprop);

      //      const double drlt = vm_.xi*q2/( pow( vm_.lambda, 2 )+vm_.xi*vm_.chi*q2 );
      weight *= pow(mB2_ / max_s_, 2. * pom_.epsilw) / prop_mx_;

      //================================================================
      // GENMXT
      //================================================================

      double dum = 0.;
      //--- vector meson mass
      switch (ifragv_) {
        case BeamMode::Elastic:
          dmxv_ = (*vm_bw_)(vm_var_);
          break;
        default:
          dmxv_ = outgoingPrimaryParticleMass(vm_var_, dum, false);
          break;
      }
      if (dmxv_ <= 0.)
        return 0.;
      //--- diffractive proton mass
      switch (ifragp_) {
        case BeamMode::Elastic:
          break;
        default:
          mY = pow(outgoingPrimaryParticleMass(difp_var_, dum, true), 2);
          break;
      }
      if (mY <= 0.)
        return 0.;

      //--- return if generated masses are bigger than CM energy
      if (mY + dmxv_ > w - 0.1)
        return 0.;

      mY2_ = mY * mY;

      //--- calculate slope parameter b
      // generate t with e**(b*t) distribution

      double b = slp_.b0 + 4. * pom_.alpha1 * log(w / slp_.wb0);
      if (ifragp_ != BeamMode::Elastic)
        b -= 4. * pom_.alpha1m * log(mY / slp_.amxb0);
      if (ifragv_ != BeamMode::Elastic)
        b -= 4. * pom_.alpha1 * log(dmxv_ / slp_.amxb0);
      b = std::max(b, 0.5);

      weight *= bmin_ / b;

      const double t = computeT(t_var_, b);
      CG_DEBUG_LOOP("DiffVM:weight") << "computed t=" << t << " GeV² for b=" << b << ".";

      //--- calculate actual minimal and maximal t for the generated masses
      // note that t here is positive!

      // Formula (E.5) from Review of Particle Properties 1992, p. III.50
      // 1: gamma, 2: p, 3: VM(+X), 4: p remnant
      // The formula for Pcm1 is altered to take the imaginary photon mass into account.

      const double inv_w = 1. / w;
      const double pcm1 = 0.5 * sqrt(pow(mB2_ + q2 - mp2_, 2) + 4. * q2 * mp2_) * inv_w;
      const double p_out = 0.5 * sqrt((mB2_ - pow(dmxv_ + mY, 2)) * (mB2_ - pow(dmxv_ - mY, 2))) * inv_w;  // pcm3
      const double t_mean =
          0.5 * ((-q2 - mp2_) * (dmxv_ * dmxv_ - mY2_) * inv_w * inv_w + mB2_ + q2 - mp2_ - dmxv_ * dmxv_ - mY2_);
      const double t_min = t_mean - 2. * pcm1 * p_out, t_max = t_mean + 2. * pcm1 * p_out;
      if (t < t_min || t > t_max)
        return 0.;

      double yhat = 0.25 * (t - t_min) / (pcm1 * p_out);
      if (yhat < 0.) {
        CG_ERROR("DiffVM:weight") << "yhat=" << yhat << " < 0.";
        yhat = 0.;
      }
      if (yhat > 1.) {
        CG_ERROR("DiffVM:weight") << "yhat=" << yhat << " > 1.";
        yhat = 1.;
      }

      //================================================================
      // GENDIF
      // /!\ in the gamma-p centre of mass frame
      //================================================================

      //--- calculate 5-vectors of diffractive states in the CMS

      Momentum p_vm_cm = p_gam_;
      p_vm_cm.lorentzBoost(p_cm_);

      // ivvm
      const double ctheta = 1. - 2. * yhat, stheta = 2. * sqrt(yhat - yhat * yhat);

      const double p_gamf = p_out * ctheta / p_vm_cm.p();
      const Momentum pt(-cos(phi_var_) * p_vm_cm.pz(),
                        sin(phi_var_) * p_vm_cm.pz(),
                        cos(phi_var_) * p_vm_cm.px() - sin(phi_var_) * p_vm_cm.py());
      const double ptf = p_out * stheta / std::hypot(p_vm_cm.pz(), pt.pz());

      p_vm_cm_ = p_gamf * p_vm_cm + ptf * pt;
      p_vm_cm_.setMass(dmxv_);

      if (sqrt(fabs(p_out * p_out - p_vm_cm_.p2())) > 0.1 * p_out)
        CG_WARNING("DiffVM:weight") << "p_out != |p_vm_cm|\n\t"
                                    << "p_out: " << p_out << ", p_vm_cm = " << p_vm_cm_ << ".";

      p_px_cm_ = -p_vm_cm_;
      p_px_cm_.setMass(mY);

      //--- calculate momentum carried by the pomeron
      // pomeron is thought to be a quasireal particle emitted by
      // the proton and absorbed by the virtual vector meson

      //std::cout << p_vmx_cm << std::endl;
      p_pom_cm_ = p_vm_cm_ - p_gam_;

      return weight;
    }

    void DiffVM::fillKinematics() {
      (*event_)[Particle::Parton1][0].setMomentum(p_gam_);
      (*event_)[Particle::OutgoingBeam1][0].setMomentum(p_gam_remn_);

      auto& op_gam = event_->getOneByRole(Particle::OutgoingBeam1);
      op_gam.setMomentum(p_gam_remn_);

      auto& pom = event_->getOneByRole(Particle::Parton2);
      Momentum p_pom_lab = p_pom_cm_;
      p_pom_lab.lorentzBoost(-p_cm_);
      (*event_)[Particle::Parton2][0].setMomentum(p_pom_lab);

      auto& op_pom = event_->getOneByRole(Particle::OutgoingBeam2);
      Momentum p_px_lab = p_px_cm_;
      p_px_lab.lorentzBoost(-p_cm_);
      (*event_)[Particle::OutgoingBeam2][0].setMomentum(p_px_lab);

      (*event_)[Particle::Intermediate][0].setMomentum(p_gam_ + p_px_lab);

      auto& vmx = event_->operator[](Particle::CentralSystem)[0];
      Momentum p_vm_lab = p_vm_cm_;
      p_vm_lab.lorentzBoost(-p_cm_);
      (*event_)[Particle::CentralSystem][0].setMomentum(p_vm_lab);
    }

    bool DiffVM::generatePhoton() {
      const Momentum p_ib = event_->oneWithRole(Particle::IncomingBeam1).momentum();
      switch (igammd_) {
        case PhotonMode::Fixed: {     // fixphot
          const double e_gamma = 3.;  // in steering card
          const double y = e_gamma / p_ib.energy();
          p_gam_ = y * p_ib;
          //p_gam_.setMass( -sqrt( fabs( p_ib.mass2()*y*y/( 1.-y ) ) ) );
          p_gam_remn_ = p_ib - p_gam_;
          p_gam_remn_.setMass(p_ib.mass());
          return true;
        } break;
        case PhotonMode::InvK: {  // genphot
          const double e_max = p_ib.p();
          const double r = exp(pho_var_ * log(min_pho_energy_ / e_max));
          if (r >= 1.)
            CG_WARNING("DiffVM:photon") << "r=" << r << " > 1.";
          p_gam_ = r * p_ib;
          p_gam_remn_ = p_ib - p_gam_;
          p_gam_remn_.setMass(p_ib.mass());
          //          std::cout << min_pho_energy_ << "/" << e_max << "->" << p_ib.mass() << std::endl;
          return true;
        } break;
        case PhotonMode::WWA:
        case PhotonMode::ABTSmith:
        case PhotonMode::AandS: {
          const auto& res = epa_calc_(pho_var_, wwa_var_);
          p_gam_ = res.pph;
          p_gam_remn_ = res.ppe;
          return res.valid;
        } break;
        default: {
          throw CG_FATAL("DiffVM:photon") << "Unsupported photon generation mode: " << igammd_ << "!";
        } break;
      }
      return false;
    }

    double DiffVM::outgoingPrimaryParticleMass(double x, double& y, bool treat) const {
      const auto& m_range = kin_.cuts().remnants.mx;
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
        CG_ERROR("DiffVM:mass") << "M=" << m << " < minimum mass=" << m_range.min() << " GeV.";
        return m_range.min();
      }
      if (m > m_range.max()) {
        CG_ERROR("DiffVM:mass") << "M=" << m << " > maximum mass=" << m_range.max() << " GeV.";
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

      return -1.;
    }

    double DiffVM::computeT(double x, double b) const {
      const auto& t_range = kin_.cuts().initial.q2();
      const double t_min = t_range.min(), t_max = t_range.max();

      double bloc = b;

      //--- generate spectrum by method of R. Lausen

      CG_DEBUG_LOOP("DiffVM:t") << "t range: " << t_range << ", b: " << b << ".";

      if (b < 0.1) {
        CG_WARNING("DiffVM:t") << "b = " << b << " < 0.1.";
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
      CG_DEBUG_LOOP("DiffVM:t") << "x=" << x << ", c0=" << c0 << ", c1=" << c1 << ", anexp=" << slp_.anexp
                                << ", bloc=" << bloc << ", t=" << t;
      return t;
    }

    std::ostream& operator<<(std::ostream& os, const DiffVM::PhotonMode& mode) {
      switch (mode) {
        case DiffVM::PhotonMode::Fixed:
          return os << "fixed energy photon";
        case DiffVM::PhotonMode::InvK:
          return os << "1/k spectrum";
        case DiffVM::PhotonMode::WWA:
          return os << "(WWA) default";
        case DiffVM::PhotonMode::ABTSmith:
          return os << "(WWA) ABT & Smith";
        case DiffVM::PhotonMode::AandS:
          return os << "(WWA) A and S";
      }
      return os;
    }

    ParametersDescription DiffVM::description() {
      auto desc = Process::description();
      desc.setDescription("Diffractive vector meson production");
      desc.add<int>("protonMode", (int)BeamMode::Elastic);
      desc.add<int>("vmMode", (int)BeamMode::Elastic);
      desc.add<int>("photonMode", (int)PhotonMode::WWA);
      desc.add<ParametersDescription>("slopeParameters", SlopeParameters::description());
      desc.add<ParametersDescription>("pomeronParameters", PomeronParameters::description());
      desc.add<ParametersDescription>("vmParameters", VectorMesonParameters::description());
      desc.add<ParametersDescription>("epaParameters", EPA::description());
      return desc;
    }

    DiffVM::SlopeParameters::SlopeParameters(const ParametersList& params) : SteeredObject(params) {
      (*this).add("b0", b0).add("wb0", wb0).add("amxb0", amxb0).add("anexp", anexp);
    }

    ParametersDescription DiffVM::SlopeParameters::description() {
      auto desc = ParametersDescription();
      desc.add<double>("b0", 4.);
      desc.add<double>("wb0", 95.);
      desc.add<double>("amxb0", 14.);
      desc.add<double>("anexp", 0.);
      return desc;
    }

    DiffVM::PomeronParameters::PomeronParameters(const ParametersList& params) : SteeredObject(params) {
      (*this).add("epsilonW", epsilw).add("epsilonM", epsilm).add("alpha1", alpha1).add("alpha1m", alpha1m);
    }

    ParametersDescription DiffVM::PomeronParameters::description() {
      auto desc = ParametersDescription();
      desc.add<double>("epsilonW", 0.225);
      desc.add<double>("epsilonM", 0.0808);
      desc.add<double>("alpha1", 0.);
      desc.add<double>("alpha1m", 0.);
      return desc;
    }

    DiffVM::VectorMesonParameters::VectorMesonParameters(const ParametersList& params) : SteeredObject(params) {
      (*this).add("lambda", lambda).add("eprop", eprop).add("xi", xi).add("chi", chi);
    }

    ParametersDescription DiffVM::VectorMesonParameters::description() {
      auto desc = ParametersDescription();
      desc.add<double>("lambda", 0.);
      desc.add<double>("eprop", 2.5);
      desc.add<double>("xi", 1.);
      desc.add<double>("chi", 1.);
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen

// register process
typedef cepgen::proc::DiffVM DiffVM;
REGISTER_PROCESS("diffvm", DiffVM);
