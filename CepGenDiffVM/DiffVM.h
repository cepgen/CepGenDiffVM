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

#include "CepGen/Physics/EPA.h"
#include "CepGen/Processes/Process.h"

namespace cepgen {
  class BreitWigner;
  namespace proc {
    /// Diffractive vector meson (photo)production as in DIFFVM \cite List:1998jz
    class DiffVM : public Process {
    public:
      explicit DiffVM(const ParametersList& = ParametersList());
      ProcessPtr clone() const override { return ProcessPtr(new DiffVM(*this)); }
      static std::string description() { return "Diffractive vector meson production"; }

      void addEventContent() override;
      void prepareKinematics() override;
      double computeWeight() override;
      void fillKinematics() override;

    private:
      bool generatePhoton();
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

      // integration variables
      double phi_var_{0.};
      double t_var_{0.};
      double pho_var_{0.};
      double wwa_var_{0.};
      double vm_var_{0.};
      double difp_var_{0.};

      /// Type of vector meson exchanged
      pdgid_t vm_pdgid_;
      /// Beam particles treatment mode
      enum class BeamMode { Elastic = 0, GluonFragmentation = -1, StandardFragmentation = 1, NucleonPionsDecay = 2 };
      BeamMode ifragp_, ifragv_;
      /// Photon generation mode
      enum class PhotonMode { Fixed = -1, InvK = 0, WWA = 1, ABTSmith = 2, AandS = 3 };
      PhotonMode igammd_;
      /// Human-readable format of a photon generation mode
      friend std::ostream& operator<<(std::ostream&, const PhotonMode&);
      struct SlopeParameters : SteeredObject<SlopeParameters> {
        explicit SlopeParameters(const ParametersList& params);

        static ParametersDescription description();

        /** \brief Slope parameter b of t distribution in GeV\f${}^{-2}\f$
           * * at CM energy \a wb0, and
           * * at mass \a amxb0 (for diffractive dissociation)
           * \note Must be positive!
           */
        double b0{0.};
        /// CM energy of \f$\gamma p\f$ system at which \f$b_0\f$ was measured, in GeV
        double wb0{0.};
        /// Mass of diffractively dissociating hadronic system for which \f$b_0\f$ was measured
        double amxb0{0.};
        /// Power law exponent
        double anexp{0.};
      } slp_;
      struct PomeronParameters : SteeredObject<PomeronParameters> {
        explicit PomeronParameters(const ParametersList&);
        /** \brief Intercept of pomeron trajectory minus 1
           * \note Controls rise of \f$\sigma_{\gamma p}\f$ with W
           */

        static ParametersDescription description();

        double epsilw{0.};
        /** \brief Intercept of pomeron trajectory minus 1
           * \note Controls \f$M_{X}\f$ spectrum
           */
        double epsilm{0.};
        /** \brief Slope alpha' of pomeron trajectory in GeV\f${}^{-2}\f$
           * \note Controls shrinkage of b slope
           */
        double alpha1{0.};
        double alpha1m{0.};
      } pom_;
      struct VectorMesonParameters : SteeredObject<VectorMesonParameters> {
        explicit VectorMesonParameters(const ParametersList& params);

        static ParametersDescription description();

        /** \brief Parameter for \f$Q^2\f$-dependence of cross section in GeV
           * \note \f$\sigma(Q^2)/\sigma(0) = 1 / \left(1 + Q^2/\Lambda^2\right)^{\rm eprop}\f$
           */
        double lambda{0.};
        /// Propagator term exponent
        double eprop{0.};
        /** \brief Parameter for \f$Q^2\f$-dependence of \f$\sigma_L/\sigma_T\f$
            * \note
            *  * \f$\frac{\sigma_L(Q^2)}{\sigma_T(Q^2)}=\frac{\xi Q^2/m^2}{1+\xi\chi Q^2/m^2}\f$ where \f$\sigma_L/\sigma_T\to \xi Q^2/m^2\f$ for low-\f$Q^2\f$, and \f$\sigma_L/\sigma_T\to 1/\chi\f$ for high-\f$Q^2\f$ ;
            *  * \f$\xi\f$ is assumed to be less than 4 (more precisely, it is assumed that \f$\sigma_L(Q^2)\f$ is always less than \f$\sigma_T(0)\f$).
            */
        double xi{0.};
        /// Purely phenomenological parameter with no theoretical justification (see \a xi)
        double chi{0.};
      } vm_;
      EPA epa_calc_;

      double bmin_{0.};
      double dmxv_{0.};
      double min_pho_energy_{0.}, max_s_{0.};
      double vm_mass_{0.}, vm_width_{0.};
      std::shared_ptr<BreitWigner> vm_bw_;
      double prop_mx_{0.};

      Momentum p_gam_, p_gam_remn_;
      Momentum p_cm_, p_pom_cm_, p_px_cm_, p_vm_cm_;
    };
  }  // namespace proc
}  // namespace cepgen

#endif
