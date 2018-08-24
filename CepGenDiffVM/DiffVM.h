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

#ifndef CepGenProcesses_DiffVM_h
#define CepGenProcesses_DiffVM_h

#include "CepGen/Processes/GenericProcess.h"

namespace CepGen {
  class BreitWigner;
  class EPA;
  namespace Process {
    class DiffVM : public GenericProcess {
    public:
      explicit DiffVM(const ParametersList&);
      ProcessPtr clone() const override { return ProcessPtr(new DiffVM(*this)); }

      void addEventContent() override;
      void setKinematics(const Kinematics&) override;
      double computeWeight() override;
      unsigned int numDimensions(const Kinematics::Mode&) const override;
      void fillKinematics() override;

    private:
      bool generatePhoton(const std::vector<double>& x);
      /// Compute the outgoing proton remnant mass
      /// \param[in] x A (collection of) random number(s) (between 0 and 1)
      /// \param[out] dw The size of the integration bin
      /// \return Mass of the outgoing proton remnant
      double outgoingPrimaryParticleMass(double x, double& y, bool treat) const;
      /// Compute the single photon virtuality for this phase space point
      /// \param[in] x Phase space point coordinate
      /// \param[in] b
      /// \return Photon virtuality, in \f${\rm GeV}^2\f$
      double computeT(double x, double b) const;

      PDG vm_pdgid_;
      enum class BeamMode {
        Elastic = 0,
        GluonFragmentation = -1,
        StandardFragmentation = 1,
        NucleonPionsDecay = 2
      } ifragp_,
          ifragv_;
      /// Photon generation mode
      enum class PhotonMode { Fixed = -1, InvK = 0, WWA = 1, ABTSmith = 2, AandS = 3 } igammd_;
      /// Human-readable format of a photon generation mode
      friend std::ostream& operator<<(std::ostream&, const PhotonMode&);
      struct SlopeParameters {
        SlopeParameters(const ParametersList& params = ParametersList());
        /// slope parameter b of t distribution in \f${\rm GeV}^{-2}\$f
        /// * at CM energy \a wb0, and
        /// * at mass \a amxb0 (for diffractive dissociation)
        /// \note Must be positive!
        double b0;
        /// CM energy of \f$\gamma p\f$ system at which \f$b_0\f$ was measured, in GeV
        double wb0;
        /// Mass of diffractively dissociating hadronic system for which \f$b_0\f$ was measured
        double amxb0;
        /// Power law exponent
        double anexp;
      } slp_;
      struct PomeronParameters {
        PomeronParameters(const ParametersList& params = ParametersList());
        /// Intercept of pomeron trajectory minus 1
        /// \note Controls rise of \f$\sigma_{\gamma p}\f$ with W
        double epsilw;
        /// Intercept of pomeron trajectory minus 1
        /// \note Controls \f$M_{X}\f$ spectrum
        double epsilm;
        /// Slope alpha' of pomeron trajectory in \f${\rm GeV}^{-2}\f$
        /// \note Controls shrinkage of b slope
        double alpha1;
        double alpha1m;
      } pom_;
      struct VectorMesonParameters {
        VectorMesonParameters(const ParametersList& params = ParametersList());
        /// Parameter for Q2-dependence of cross section in GeV
        /// \note \f$\sigma(Q^2)/\sigma(0) = 1 / \left(1 + Q^2/\Lambda^2\right)^{\rm eprop}\f$
        double lambda;
        /// Propagator term exponent
        double eprop;
        /// Parameter for \f$Q^2\f$-dependence of \f$\sigma_L/\sigma_T\f$
        /// \note
        ///  * \f$\frac{\sigma_L(Q^2)}{\sigma_T(Q^2)}=\frac{\xi Q^2/m^2}{1+\xi\chi Q^2/m^2}\f$ where \f$\sigma_L/\sigma_T\to \xi Q^2/m^2\f$ for low-\f$Q^2\f$, and \f$\sigma_L/\sigma_T\to 1/\chi\f$ for high-\f$Q^2\f$ ;
        ///  * \f$\xi\f$ is assumed to be less than 4 (more precisely, it is assumed that \f$\sigma_L(Q^2)\f$ is always less than \f$\sigma_T(0)\f$).
        double xi;
        /// Purely phenomenological parameter with no theoretical justification (see \a xi)
        double chi;
      } vm_;

      double bmin_;
      double dmxv_;
      double min_pho_energy_, max_s_;
      double vm_mass_, vm_width_;
      std::shared_ptr<BreitWigner> vm_bw_;
      std::shared_ptr<EPA> epa_calc_;
      double prop_mx_;
      unsigned short ndim_;

      Particle::Momentum p_gam_, p_gam_remn_;
      Particle::Momentum p_cm_, p_pom_cm_, p_px_cm_, p_vm_cm_;
    };
  }  // namespace Process
}  // namespace CepGen

#endif
