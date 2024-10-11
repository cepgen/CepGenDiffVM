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

#include <cmath>

#include "CepGenDiffVM/BreitWigner.h"

namespace cepgen::utils {
  BreitWigner::BreitWigner(double mean, double gamma, const Limits& energy_range)
      : mean_(mean), gamma_(gamma), energy_range_(energy_range) {}

  double BreitWigner::operator()(double x) const {
    if (const double val = mean_ + 0.5 * gamma_ * std::tan((2. * x - 1.) * M_PI_2); energy_range_.contains(val))
      return val;
    return -1.;
  }
}  // namespace cepgen::utils
