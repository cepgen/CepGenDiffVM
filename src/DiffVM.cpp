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
    DiffVM::DiffVM() : GenericProcess("diffvm", "Exclusive vector meson production"), bmin_(0.), dmxv_(0.), dmxp_(0.) {
      if (ifragp_ == ProtonMode::Elastic) {
        bmin_ = slp_.b0 + 4. * pom_.alpha1 * log(cuts_.cuts.central.mass_single.min() / slp_.wb0);
        dmxp_ = mp_;
      } else {
        if (ifragv_ == VectorMode::Elastic)
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
      const double yhat = 0, pout = 0.;

      //===== GENDIF

      const double ctheta = 1. - 2. * yhat, stheta = 2. * sqrt(yhat - yhat * yhat);

      //--- calculate 5-vectors of diffractive states in the CMS

      const double p_gamf = pout * ctheta / p_vm_.p();
      const double phi = 2. * M_PI * x(0);
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

    unsigned int DiffVM::numDimensions(const Kinematics::Mode&) const { return 1; }

    void DiffVM::fillKinematics() {}

    double DiffVM::computeOutgoingPrimaryParticlesMasses(double x, double& y) {
      const auto& m_range = cuts_.cuts.remnants.mass_single;
      const double mmin = m_range.min(), mmax = m_range.max();
      double lmin = 0., delta = 0., m2min = 0., fact = 0., m2 = 0.;
      if (fabs(pom_.epsilm) < 1.e-3) {
        lmin = 2. * log(mmin);
        delta = 2. * log(mmax / mmin);
      } else {
        m2min = pow(mmin, -2 * pom_.epsilm);
        fact = pow(mmax, -2 * pom_.epsilm) - m2min;
      }
      unsigned short i = 0;
      while (true) {
        if (fabs(pom_.epsilm) < 1.e-3)
          //--- basic spectrum: 1/M^2
          m2 = exp(x * delta + lmin);
        else
          //--- basic spectrum: 1/M^2(1+epsilon)
          m2 = pow(fact * x + m2min, -1. / pom_.epsilm);

        if (m2 < mmin * mmin) {
          CG_ERROR("DiffVM:mass") << "M2 =" << m2 << " < MMIN**2 = " << mmin * mmin << ".\n\t"
                                  << " LMIN = " << lmin << ", DELTA = " << delta << ",\n\t"
                                  << " M2MIN = " << m2min << ", FACT = " << fact << ",\n\t"
                                  << " EPSILM = " << pom_.epsilm << ".";
          m2 = mmin * mmin;
        } else if (m2 > mmax * mmax) {
          CG_ERROR("DiffVM:mass") << "M2 =" << m2 << " > MMIN**2 = " << mmin * mmin << ".";
          m2 = mmax * mmax;
        }

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

        if (1.60 * rand() < y || ++i >= 100)
          break;

        /*//--- new version: 1/M_x^2 spectrum
        if ( m2 >= 2. ) y = 1.;
        else            y = 1.-1.1815*pow( m2-2., 2 );
        if ( rand() < y || ++i >= 100 )
          break;*/
      }

      if (y > 1.6)
        CG_WARNING("DiffVM:mass") << "Y = " << y << " for M2 = " << m2 << ".";

      CG_DEBUG_LOOP("DiffVM:mass") << "Number of iterations: " << i << ".";
      if (i > 100)
        CG_WARNING("DiffVM:mass") << "More than 100 iterations!\n\t"
                                  << "MMIN: " << mmin << " MMAX: " << mmax << " M2 = " << m2 << ".";

      return sqrt(m2);
    }
  }  // namespace Process
}  // namespace CepGen
