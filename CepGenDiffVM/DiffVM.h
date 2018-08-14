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
      void fillKinematics(bool) override;
    };
  }  // namespace Process
}  // namespace CepGen

#endif
