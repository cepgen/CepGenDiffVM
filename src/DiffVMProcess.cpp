/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2024  Laurent Forthomme
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

#include <CepGen/Core/Exception.h>
#include <CepGen/Event/Event.h>
#include <CepGen/Modules/ProcessFactory.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Physics/ParticleProperties.h>
#include <CepGen/Process/Process.h>

#include <cmath>

#include "CepGenDiffVM/BreitWigner.h"
#include "CepGenDiffVM/EPA.h"

using namespace cepgen;

/// Diffractive vector meson (photo)production as in DIFFVM \cite List:1998jz
class DiffVMProcess : public cepgen::proc::Process {
public:
  explicit DiffVMProcess(const ParametersList& params)
      : proc::Process(params),
        vm_pdgid_(steer<ParticleProperties>("vmFlavour").pdgid),
        ifragp_(steerAs<int, BeamMode>("protonMode")),
        ifragv_(steerAs<int, BeamMode>("vmMode")),
        igammd_(steerAs<int, PhotonMode>("photonMode")),
        slope_parameters_(steer<ParametersList>("slopeParameters")),
        pomeron_parameters_(steer<ParametersList>("pomeronParameters")),
        vm_parameters_(steer<ParametersList>("vmParameters")),
        epa_calculator_(steer<ParametersList>("epaParameters")) {}

  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new DiffVMProcess(*this)); }

  static ParametersDescription description() {
    auto desc = proc::Process::description();
    desc.setDescription("Diffractive vector meson production");
    desc.addAs<int, BeamMode>("protonMode", BeamMode::Elastic);
    desc.addAs<int, BeamMode>("vmMode", BeamMode::Elastic);
    desc.addAs<int, PhotonMode>("photonMode", PhotonMode::WWA);
    desc.add<double>("fixphot", 3.);
    desc.add<ParametersDescription>("slopeParameters", SlopeParameters::description());
    desc.add<ParametersDescription>("pomeronParameters", PomeronParameters::description());
    desc.add<ParametersDescription>("vmParameters", VectorMesonParameters::description());
    desc.add<ParametersDescription>("epaParameters", EPA::description());
    return desc;
  }

  void addEventContent() override {
    Process::setEventContent({{Particle::Role::IncomingBeam1, {PDG::electron}},
                              {Particle::Role::IncomingBeam2, {PDG::proton}},
                              {Particle::Role::Parton1, {PDG::photon}},
                              {Particle::Role::Parton2, {PDG::pomeron}},
                              {Particle::Role::OutgoingBeam1, {PDG::electron}},
                              {Particle::Role::OutgoingBeam2, {PDG::proton}},
                              {Particle::Role::CentralSystem, {vm_pdgid_}}});
  }
  void prepareKinematics() override {
    mY2() = pB().mass2();

    //const auto& w_limits = kinematics().cuts().central.mass_single;//FIXME
    const auto& w_limits = Limits{0., 50.};

    //--- variables mapping
    defineVariable(m_phi_var_, Mapping::linear, {0., 2. * M_PI}, "azumuthal_angle");  //TODO factor 2*pi? to be checked
    defineVariable(m_t_spectrum_, Mapping::linear, {0., 1.}, "t_spectrum_variable");
    defineVariable(m_photon_rel_fraction_, Mapping::linear, {0., 1.}, "photon_fractional_momentum");
    defineVariable(m_diffractive_mass_, Mapping::linear, {0., 1.}, "diffractive_proton_mass");
    if (igammd_ >= PhotonMode::WWA) {  // use the equivalent photon approximation for the photon part
      epa_calculator_.init(pA(), pB(), kinematics().cuts().initial.q2.at(0), w_limits);
      defineVariable(m_wwa_var_, Mapping::linear, {0., 1.}, "epa_variable");
    }

    const auto compute_bmin = [this](double mass_min) {
      return std::max(
          0.5, slope_parameters_.b0 + 4. * pomeron_parameters_.alpha1 * std::log(mass_min / slope_parameters_.wb0));
    };
    if (ifragp_ == BeamMode::Elastic) {
      if (!w_limits.hasMin())
        throw CG_FATAL("DiffVM") << "You must specify a lower limit to W(gamma,p)!\n\t"
                                 << "Current limits: " << w_limits << ".";
      bmin_ = compute_bmin(w_limits.min());
    } else {
      if (ifragv_ == BeamMode::Elastic)
        bmin_ = compute_bmin(slope_parameters_.amxb0);
      else
        bmin_ = compute_bmin(4. * std::pow(slope_parameters_.amxb0, 2) * inverseSqrtS());
    }
    if (ifragv_ == BeamMode::Elastic)
      defineVariable(m_vm_var_, Mapping::linear, {0., 1.}, "VMvar");
    CG_DEBUG("DiffVMProcess:prepareKinematics") << "Minimum b slope: " << bmin_ << ".";

    min_pho_energy_ = 0.25 * std::pow(w_limits.min(), 2) / pB().p();
    max_s_ = std::pow(w_limits.max(), 2);

    if (vm_parameters_.lambda2 <= 0.)
      vm_parameters_.lambda2 = event().oneWithRole(Particle::Role::CentralSystem).momentum().mass2();

    prop_mx_ = vm_parameters_.computeMX(kinematics().cuts().initial.q2.at(0).min());

    const Particle& vm = event()(Particle::Role::CentralSystem)[0];
    vm_mass_ = vm.momentum().mass(), vm_width_ = PDG::get().width(vm.pdgId());

    //--- mass range for VM generation
    double min_vm_mass = -1., max_vm_mass = -1.;
    const auto& invm_limits = kinematics().cuts().central.mass_sum;
    if (invm_limits.valid()) {
      min_vm_mass = invm_limits.min();
      max_vm_mass = invm_limits.max();
    } else {
      min_vm_mass = vm_mass_ - 3. * vm_width_;
      max_vm_mass = vm_mass_ + 10. * vm_width_;
      if (vm.pdgId() == 100113 /*rho1450_0*/ || vm.pdgId() == 30113 /*rho1700_0*/)
        min_vm_mass = std::max(min_vm_mass, 1.2);
      else if (vm.pdgId() == 10333 /*h1380_1*/)
        min_vm_mass = std::max(min_vm_mass, 1.4);
    }
    vm_bw_.reset(new utils::BreitWigner(vm_mass_, vm_width_, {min_vm_mass, max_vm_mass}));
  }
  double computeWeight() override;
  void fillKinematics() override {
    q2().lorentzBoost(-p_cm_);
    pY().lorentzBoost(-p_cm_);
    event().oneWithRole(Particle::Role::Intermediate).setMomentum(q1() + q2());
    event().oneWithRole(Particle::Role::CentralSystem).setMomentum(Momentum(p_vm_cm_).lorentzBoost(-p_cm_));
  }

private:
  spdgid_t vm_pdgid_;  ///< Type of vector meson exchanged
  /// Beam particles treatment mode
  enum class BeamMode { Elastic = 0, GluonFragmentation = -1, StandardFragmentation = 1, NucleonPionsDecay = 2 };
  BeamMode ifragp_, ifragv_;
  enum class PhotonMode { Fixed = -1, InvK = 0, WWA = 1, ABTSmith = 2, AandS = 3 } igammd_;  ///< Photon generation mode
  /// Human-readable format of a photon generation mode
  friend std::ostream& operator<<(std::ostream& os, const PhotonMode& pm) {
    switch (pm) {
      case PhotonMode::Fixed:
        return os << "fixed energy photon";
      case PhotonMode::InvK:
        return os << "1/k spectrum";
      case PhotonMode::WWA:
        return os << "(WWA) default";
      case PhotonMode::ABTSmith:
        return os << "(WWA) ABT & Smith";
      case PhotonMode::AandS:
        return os << "(WWA) A and S";
    }
    return os;
  }
  struct SlopeParameters : SteeredObject<SlopeParameters> {
    explicit SlopeParameters(const ParametersList& params) : SteeredObject(params) {
      (*this).add("b0", b0).add("wb0", wb0).add("amxb0", amxb0).add("anexp", anexp);
    }
    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.add<double>("b0", 4.);
      desc.add<double>("wb0", 95.);
      desc.add<double>("amxb0", 14.);
      desc.add<double>("anexp", 0.);
      return desc;
    }
    /// Slope parameter b of t distribution in GeV\f${}^{-2}\f$
    ///  * at CM energy \a wb0, and
    ///  * at mass \a amxb0 (for diffractive dissociation)
    /// \note Must be positive!
    double b0{0.};
    double wb0{0.};    ///< CM energy of \f$\gamma p\f$ system at which \f$b_0\f$ was measured, in GeV
    double amxb0{0.};  ///< Mass of diffractively dissociating hadronic system for which \f$b_0\f$ was measured
    double anexp{0.};  ///< Power law exponent
  } slope_parameters_;

  struct PomeronParameters : SteeredObject<PomeronParameters> {
    explicit PomeronParameters(const ParametersList& params) : SteeredObject(params) {
      (*this).add("epsilonW", epsilw).add("epsilonM", epsilm).add("alpha1", alpha1).add("alpha1m", alpha1m);
    }
    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.add<double>("epsilonW", 0.225);
      desc.add<double>("epsilonM", 0.0808);
      desc.add<double>("alpha1", 0.);
      desc.add<double>("alpha1m", 0.);
      return desc;
    }
    /// Intercept of pomeron trajectory minus 1
    /// \note Controls rise of \f$\sigma_{\gamma p}\f$ with W
    double epsilw{0.};
    /// Intercept of pomeron trajectory minus 1
    /// \note Controls \f$M_{X}\f$ spectrum
    double epsilm{0.};
    /// Slope alpha' of pomeron trajectory in GeV\f${}^{-2}\f$
    /// \note Controls shrinkage of b slope
    double alpha1{0.};
    double alpha1m{0.};
  } pomeron_parameters_;

  struct VectorMesonParameters : SteeredObject<VectorMesonParameters> {
    explicit VectorMesonParameters(const ParametersList& params)
        : SteeredObject(params),
          lambda2(std::pow(steer<double>("lambda"), 2)),
          eprop(steer<double>("eprop")),
          xi(steer<double>("xi")),
          chi(steer<double>("chi")) {}

    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.add<double>("lambda", 0.).setDescription("Q^2-dependence parameter of cross section in GeV");
      desc.add<double>("eprop", 2.5).setDescription("Propagator term exponent");
      desc.add<double>("xi", 1.).setDescription("Q^2-dependence parameter of sigma_L/sigma_T");
      desc.add<double>("chi", 1.).setDescription("Purely phenomenological parameter with no theoretical justification");
      return desc;
    }
    double computeMX(double q2_min) const {
      return std::max(1., xi * q2_min / (lambda2 + xi * chi * q2_min) / std::pow(1. + q2_min / lambda2, eprop));
    }
    /// Parameter for \f$Q^2\f$-dependence of cross section in GeV^2
    /// \note \f$\sigma(Q^2)/\sigma(0) = 1 / \left(1 + Q^2/\Lambda^2\right)^{\rm eprop}\f$
    double lambda2{0.};
    double eprop{0.};  ///< Propagator term exponent
    /// Parameter for \f$Q^2\f$-dependence of \f$\sigma_L/\sigma_T\f$
    /// \note
    ///  * \f$\frac{\sigma_L(Q^2)}{\sigma_T(Q^2)}=\frac{\xi Q^2/m^2}{1+\xi\chi Q^2/m^2}\f$ where \f$\sigma_L/\sigma_T\to \xi Q^2/m^2\f$ for low-\f$Q^2\f$, and \f$\sigma_L/\sigma_T\to 1/\chi\f$ for high-\f$Q^2\f$ ;
    ///  * \f$\xi\f$ is assumed to be less than 4 (more precisely, it is assumed that \f$\sigma_L(Q^2)\f$ is always less than \f$\sigma_T(0)\f$).
    double xi{0.};
    double chi{0.};  ///< Purely phenomenological parameter with no theoretical justification (see \a xi)
  } vm_parameters_;
  EPA epa_calculator_;

  double bmin_{0.};
  double dmxv_{0.};
  double min_pho_energy_{0.}, max_s_{0.};
  double vm_mass_{0.}, vm_width_{0.};
  std::shared_ptr<utils::BreitWigner> vm_bw_;
  double prop_mx_{0.};

  Momentum p_cm_, p_vm_cm_;

  // integration variables
  double m_phi_var_{0.};
  double m_t_spectrum_{0.};
  double m_photon_rel_fraction_{0.};
  double m_wwa_var_{0.};
  double m_vm_var_{0.};
  double m_diffractive_mass_{0.};
};

double DiffVMProcess::computeWeight() {
  const auto generate_photon = [this]() {  // GENGAM in DiffVM
    switch (igammd_) {
      case PhotonMode::Fixed: {
        const double e_gamma = steer<double>("fixphot");
        const double y = e_gamma / pA().energy();
        q1() = y * pA();
        //q1().setMass(-std::sqrt(fabs(pA().mass2() * y * y / (1. - y))));
        pX() = (pA() - q1()).setMass(pA().mass());
        return true;
      } break;
      case PhotonMode::InvK: {  // genphot
        const double e_max = pA().p();
        const double r = std::exp(m_photon_rel_fraction_ * std::log(min_pho_energy_ / e_max));
        if (r >= 1.)
          CG_WARNING("DiffVMProcess:photon") << "r=" << r << " > 1.";
        q1() = r * pA();
        pX() = (pA() - q1()).setMass(pA().mass());
        return true;
      } break;
      case PhotonMode::WWA:
      case PhotonMode::ABTSmith:
      case PhotonMode::AandS: {
        const auto& res = epa_calculator_(m_photon_rel_fraction_, m_wwa_var_);
        q1() = res.pph;
        pX() = res.ppe;
        return res.valid;
      } break;
      default: {
        throw CG_FATAL("DiffVMProcess:photon") << "Unsupported photon generation mode: " << igammd_ << "!";
      } break;
    }
    return false;
  };
  if (!generate_photon())
    return 0.;

  const double q2_value = q1().mass2();
  if (!kinematics().cuts().initial.q2.at(0).contains(q2_value))
    return 0.;

  //--- determine gamma*-p energy
  p_cm_ = q1() + pB();
  const auto mb2 = p_cm_.energy2(), mb = p_cm_.energy();

  double weight = 1.;
  weight /= std::pow(1. + q2_value / vm_parameters_.lambda2, vm_parameters_.eprop);  // weight of the virtual VM
  //const double drlt =
  //    vm_parameters_.xi * q2_value / (vm_parameters_.lambda2 + vm_parameters_.xi * vm_parameters_.chi * q2_value);
  weight *= std::pow(mb2 / max_s_, 2. * pomeron_parameters_.epsilw) / prop_mx_;

  const auto generate_outgoing_particle_mass = [this](double x) {  // GENMXT in DiffVM
    const auto& mass_range = kinematics().cuts().remnants.mx;
    if (std::fabs(pomeron_parameters_.epsilm) < 1.e-3)  // basic spectrum: 1/M^2
      return mass_range.trim(mass_range.min() * std::pow(mass_range.max() / mass_range.min(), x));
    // basic spectrum: 1/M^2(1+epsilon)
    const double m2min = std::pow(mass_range.min(), -2. * pomeron_parameters_.epsilm);
    const double fact = std::pow(mass_range.max(), -2. * pomeron_parameters_.epsilm) - m2min;
    return mass_range.trim(std::sqrt(std::pow(fact * x + m2min, -1. / pomeron_parameters_.epsilm)));
  };
  const auto compute_vm_mass = [this, &generate_outgoing_particle_mass](double x) {
    switch (ifragv_) {
      case BeamMode::Elastic:
        return (*vm_bw_)(x);
      default:
        return generate_outgoing_particle_mass(x);
    }
  };
  const auto compute_diffractive_mass = [this, &generate_outgoing_particle_mass](double x) {
    switch (ifragp_) {
      case BeamMode::Elastic:
        return mB();
      default: {
        const auto mass = generate_outgoing_particle_mass(x);
        const auto generate_y = [](double mass2) -> double {  // old version with enhancements in lower mass region
          if (mass2 >= 4.)
            return 1.;
          if (mass2 >= 3.1)
            return 1.64 - 0.16 * mass2;
          if (mass2 >= 2.65)
            return mass2 * (0.47 - 0.42 * std::pow(mass2 - 2.65, 2));
          if (mass2 >= 2.25)
            return mass2 * (0.47 + 0.46 * std::pow(mass2 - 2.65, 2));
          if (mass2 >= 2.02)
            return mass2 * (0.76 - 2.69 * std::pow(mass2 - 2.02, 2));
          if (mass2 >= 1.72)
            return mass2 * (0.76 - 1.98 * std::pow(mass2 - 2.02, 2));
          return 1.05 * (mass2 - 1.165);
        };
        const auto generate_y_new = [](double mass2) -> double {  // new version: 1/M_x^2 spectrum
          return mass2 < 2. ? 1. - 1.1815 * std::pow(mass2 - 2., 2) : 1.;
        };
        if (const auto dx = generate_y(mass * mass); 1.60 * rand() >= dx)
          return -1.;
        return mass;
      }
    }
  };
  if (dmxv_ = compute_vm_mass(m_vm_var_); dmxv_ <= 0.)  // vector meson mass
    return 0.;
  double mY = 0;
  if (mY = compute_diffractive_mass(m_diffractive_mass_); mY <= 0.)  // diffractive proton mass
    return 0.;
  if (mY + dmxv_ > mb - 0.1)  // skip if generated masses are bigger than CM energy
    return 0.;

  mY2() = mY * mY;

  //--- calculate slope parameter b
  // generate t with e**(b*t) distribution

  double b = slope_parameters_.b0 + 4. * pomeron_parameters_.alpha1 * std::log(mb / slope_parameters_.wb0);
  if (ifragp_ != BeamMode::Elastic)
    b -= 4. * pomeron_parameters_.alpha1m * std::log(mY / slope_parameters_.amxb0);
  if (ifragv_ != BeamMode::Elastic)
    b -= 4. * pomeron_parameters_.alpha1 * std::log(dmxv_ / slope_parameters_.amxb0);
  b = std::max(b, 0.5);

  weight *= bmin_ / b;

  const auto compute_t = [this](double x, double b) {  /// single photon virtuality for this phase space point
    const auto& t_range = kinematics().cuts().initial.q2.at(0);
    if (!t_range.valid()) {
      CG_ERROR("DiffVMProcess:t") << "t range: " << t_range << "=> return " << t_range.min() << ".";
      return t_range.min();
    }
    {  // ensure b is within a correct range
      if (b < 0.1)
        CG_WARNING("DiffVMProcess:t") << "b = " << b << " < 0.1.";
      b = std::max(b, 0.1);
    }
    CG_DEBUG_LOOP("DiffVMProcess:t") << "t range: " << t_range << ", b: " << b << ".";
    // generate spectrum by method of R. Lausen
    if (slope_parameters_.anexp < 1.) {  // power law exponent is 0 or illegal => generate pure exp(bt) spectrum
      if (b * t_range.range() >= 25.)    // method 1
        return t_range.min() - std::log(x) / b;
      return t_range.min() - std::log(1. - x * (1. - std::exp(b * t_range.range()))) / b;  // method 2
    }
    //--- new 16.5.07 BL:
    // Generate mixed exp(bt)/power law spectrum
    // d sigma/d t = exp (-n*ln(-bt/n+1)) = (-bt/n+1)^-n
    // Limit for small bt: exp (bt + c t^2) with c=b^2/2n
    // Limit for large bt>>n: t^-n
    const auto expo = 1. - slope_parameters_.anexp;
    const auto c_range =
        (t_range * b + slope_parameters_.anexp).compute([&expo](double ext) -> double { return std::pow(ext, expo); });
    const auto t = -(slope_parameters_.anexp - std::pow(c_range.x(x), 1. / (1. - slope_parameters_.anexp))) / b;
    CG_DEBUG_LOOP("DiffVMProcess:t") << "x=" << x << ", c range=" << c_range << ", anexp=" << slope_parameters_.anexp
                                     << ", b=" << b << " -> t=" << t;
    return t;
  };
  const double t = compute_t(m_t_spectrum_, b);
  CG_DEBUG_LOOP("DiffVMProcess:weight") << "computed t=" << t << " GeV² for b=" << b << ".";

  //--- calculate actual minimal and maximal t for the generated masses
  // note that t here is positive!

  // Formula (E.5) from Review of Particle Properties 1992, p. III.50
  // 1: gamma, 2: p, 3: VM(+X), 4: p remnant
  // The formula for Pcm1 is altered to take the imaginary photon mass into account.

  const double inv_w = 1. / mb;
  const double pcm1 = 0.5 * std::sqrt(std::pow(mb2 + q2_value - mp2_, 2) + 4. * q2_value * mp2_) * inv_w;
  const double p_out =
      0.5 * std::sqrt((mb2 - std::pow(dmxv_ + mY, 2)) * (mb2 - std::pow(dmxv_ - mY, 2))) * inv_w;  // pcm3
  const double t_mean = 0.5 * ((-q2_value - mp2_) * (dmxv_ * dmxv_ - mY2()) * inv_w * inv_w + mb2 + q2_value - mp2_ -
                               dmxv_ * dmxv_ - mY2());
  const auto t_range = Limits{t_mean - 2. * pcm1 * p_out, t_mean + 2. * pcm1 * p_out};
  if (!t_range.contains(t))
    return 0.;

  const auto yhat = Limits{0., 1.}.trim(0.25 * (t - t_range.min()) / (pcm1 * p_out));

  //================================================================
  // GENDIF
  // /!\ in the gamma-p centre of mass frame
  //================================================================

  //--- calculate 5-vectors of diffractive states in the CMS
  auto p_vm_cm = Momentum(q1()).lorentzBoost(p_cm_);

  const auto cos_theta = 1. - 2. * yhat, sin_theta = 2. * std::sqrt(yhat - yhat * yhat);  // ivvm
  const double p_gamf = p_out * cos_theta / p_vm_cm.p();
  const auto pt = Momentum(-std::cos(m_phi_var_) * p_vm_cm.pz(),
                           +std::sin(m_phi_var_) * p_vm_cm.pz(),
                           +std::cos(m_phi_var_) * p_vm_cm.px() - std::sin(m_phi_var_) * p_vm_cm.py());
  const double ptf = p_out * sin_theta / std::hypot(p_vm_cm.pz(), pt.pz());

  p_vm_cm_ = (p_gamf * p_vm_cm + ptf * pt).setMass(dmxv_);

  /*if (std::hypot(p_out, p_vm_cm_.p()) > 0.1 * p_out)
    CG_WARNING("DiffVMProcess:weight") << "p_out != |p_vm_cm|\n\t"
                                       << "p_out: " << p_out << ", p_vm_cm = " << p_vm_cm_ << ".";*/

  pY() = (-p_vm_cm_).setMass(mY);
  q2() = (p_vm_cm_ - q1());  // mom. carried by the pomeron (quasireal proton-emitted particle absorbed by virtual VM)
  CG_LOG << event();
  return weight;
}

REGISTER_PROCESS("diffvm", DiffVMProcess);
