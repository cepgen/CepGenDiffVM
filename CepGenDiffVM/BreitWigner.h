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

#ifndef CepGenDiffVM_BreitWigner_h
#define CepGenDiffVM_BreitWigner_h

#include <CepGen/Utils/Limits.h>

namespace cepgen::utils {
  /// A Breit-Wigner/Cauchy distribution generator
  class BreitWigner {
  public:
    explicit BreitWigner(double mean = 0., double gamma = 0., const Limits& energy_range = {});

    const Limits energyRange() const { return energy_range_; }

    double operator()(double x) const;  ///< Shoot a value according to parameterisation

  private:
    const double mean_;          ///< Mean of distribution
    const double gamma_;         ///< Width of distribution
    const Limits energy_range_;  ///< Allowed energy range
  };
}  // namespace cepgen::utils

#endif
