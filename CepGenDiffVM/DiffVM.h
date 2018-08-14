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
  namespace Process {
    class DiffVM : public GenericProcess {
    public:
      explicit DiffVM();
      ProcessPtr clone() const override { return ProcessPtr(new DiffVM(*this)); }

      void addEventContent() override;
      double computeWeight() override;
      unsigned int numDimensions(const Kinematics::Mode&) const override;
      void fillKinematics() override;

    private:
      /// Compute the ougoing proton remnant mass
      /// \param[in] x A random number (between 0 and 1)
      /// \param[out] dw The size of the integration bin
      /// \return Mass of the outgoing proton remnant
      double computeOutgoingPrimaryParticlesMasses(double x, double& y);

      enum class ProtonMode { Elastic = 0, GluonFragmentation = -1, StandardFragmentation = 1, NucleonPionsDecay = 2 };
      ProtonMode ifragp_;
      struct SlopeParameters {
        SlopeParameters() : b0(4.), wb0(95.), amxb0(14.), anexp(0.) {}
        /// slope parameter b of t distribution in GeV^-2
        /// * at CM energy wb0, and
        /// * at mass amxb0 (for diffractive dissociation)
        /// \note Must be positive!
        double b0;
        /// CM energy of gamma-p system at which b0 was measured, in GeV
        double wb0;
        /// Mass of diffractively dissociating hadronic system for which b0 was measured
        double amxb0;
        /// Power law exponent
        double anexp;
      };
      SlopeParameters slp_;
      struct PomeronParameters {
        PomeronParameters() : epsilw(0.225), epsilm(0.0808), alpha1(0.) {}
        /// Intercept of pomeron trajectory minus 1
        /// \note Controls rise of sigma_gammap with W
        double epsilw;
        /// Intercept of pomeron trajectory minus 1
        /// \note Controls Mx spectrum
        double epsilm;
        /// Slope alpha' of pomeron trajectory in GeV^-2
        /// \note Controls shrinkage of b slope
        double alpha1;
      };
      PomeronParameters pom_;
      enum class VectorMode { Elastic = 0, StandardFragmentation = 1, VectorPionsDecay = 2 };
      VectorMode ifragv_;

      double bmin_;
      Particle::Momentum p_vm_;
      double dmxv_, dmxp_;
    };
  }  // namespace Process
}  // namespace CepGen

#endif
