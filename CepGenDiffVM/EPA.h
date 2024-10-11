/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2024  Laurent Forthomme
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

#ifndef CepGenDiffVM_EPA_h
#define CepGenDiffVM_EPA_h

#include <CepGen/Core/SteeredObject.h>
#include <CepGen/Physics/Momentum.h>
#include <CepGen/Utils/Limits.h>

#include <array>
#include <iosfwd>

namespace cepgen {
  /// \brief An Equivalent Photon Approximation calculator
  /// \author Benno List
  /// \author Thomas Jansen
  /// \date 1993-1995
  class EPA : public SteeredObject<EPA> {
  public:
    /// Default constructor
    explicit EPA(const ParametersList&);

    static ParametersDescription description();

    /// Output format for a given approximation
    struct Result {
      bool valid{false};
      Momentum pph;   ///< 5-vector of photon in lab system (5th component is \f${\rm sign}(q^2)\cdot\sqrt{|q^2|}\f$)
      Momentum ppe;   ///< 5-vector of scattered electron
      double q2{0.};  ///< Virtuality of photon (positive!): \f$Q^2 = -q^2\f$
      /// Photon helicity:
      /// - 0: longitudinal,
      /// - \f$\pm\f$1: transverse polarization
      short heli{0};
    };

    /// Photon generation mode
    enum struct Mode {
      wwa = 1,                      ///< WWA approximation (including e-mass effect and longitudinal flux). Recommended
      transverse = 2,               ///< transverse spectrum
      transverse_longitudinal = 3,  ///< transverse & longitudinal spectrum
      transverse_longitudinal_pframe = 4  ///< as @a transverse_longitudinal, but flux in proton rest frame
    };
    friend std::ostream& operator<<(std::ostream& os, const Mode& mode);

    /// Initialize histograms, constants, kinematical bounds
    void init(const Momentum& pel, const Momentum& ppr, const Limits& q2_range, const Limits& w_range);

    /// Generate one event with unweighted photon and electron
    /** \note
       * 1. according to WWA:
       *    - transversal photonspectrum (\f$Q^2 \to 0\f$):
       *      \f$ P(y, Q^2) = \frac{\alpha}{2\pi}\frac{1}{Q^2 y}\left(2(1-y)\left(1-\frac{Q^2_{\rm min}}{Q^2}\right)+y^2\right),\f$
       *    - longitudinal photonspectrum (\f$Q^2 \to 0\f$):
       *      \f$ P(y, Q^2) = \frac{\alpha}{2\pi}\frac{1}{Q^2 y}2(1-y).\f$
       * 2. full transversal photonspectrum given by: \cite Abt:epa,\cite Smith:1993jqa,\cite Smith:1994ef
       * 3. full transversal and longitudinal spectrum by \cite Abt:epa
       *    - calculate integrated factor over the spectrum:
       *       kinematical bounds \f$Y_{\rm min}\f$,  \f$Y_{\rm max} (W_{\rm min})\f$, \f$Q^2_{\rm min}\f$, \f$Q^2_{\rm max} (Q^2_{\rm cutoff})\f$.
       */
    Result operator()(double x1, double x2) const;

  private:
    /// Generate the helicity of a photon
    /// \author Benno List
    /// \date 27 May 1993
    /// \param[in] longfr Fraction of longitudinally polarized photons
    /// \note IHELI in DIFFVM
    /// \return Helicity of the photon (\f$\pm\f$1: transverse photon, 0: longitudinal photon)
    short helicity(double longfr) const;
    /// alpha/2pi
    static const double ALPHARED;
    const Mode mode_;
    const Limits y_range_;
    Limits dy_range_;
    Limits q2_range_, w_range_;
    /// 5-vector of beam electron
    Momentum pel_;
    /// 5-vector of beam proton
    Momentum ppr_;
    double s_{0.}, w12_{0.}, elpr_{0.}, eel_{0.};
    mutable std::array<unsigned int, 2> num_errors_{0, 0};
    const std::array<unsigned int, 2> max_errors_{10, 500};
    mutable double epa_max_{-1.};
  };
}  // namespace cepgen

#endif
