#ifndef CepGenDiffVM_EPA_h
#define CepGenDiffVM_EPA_h

#include <ostream>

#include "CepGen/Physics/Momentum.h"
#include "CepGen/Utils/Limits.h"

namespace CepGen {
  /// \author Benno List
  /// \author Thomas Jansen
  /// \date 1993-1995
  class EPA {
  public:
    /// Photon generation mode
    enum struct Mode {
      wwa = 1,                      ///< WWA approximation (including e-mass effect and longitudinal flux). Recommended
      transverse = 2,               ///< transverse spectrum
      transverse_longitudinal = 3,  ///< transverse & longitudinal spectrum
      transverse_longitudinal_pframe = 4  ///< as 3, but flux in p rest frame
    };
    friend std::ostream& operator<<(std::ostream& os, const Mode& mode);

    EPA(const Mode& mode = Mode::wwa);
    struct Result {
      Result() : valid(false), q2(0.), heli(0) {}
      bool valid;
      /// 5-vector of photon in H1 lab. system (5th component is sign (q2)*sqrt (abs (q2)))
      Momentum pph;
      /// 5-vector of scattered electron
      Momentum ppe;
      /// Virtuality of photon (positive!): Q2 = -q2
      double q2;
      /// Photon helicity:
      /// - 0: longitudinal,
      /// - +/-1: transverse polarization
      short heli;
    };

    /// Initialize histograms, constants, kinematical bounds
    void init(const Momentum& pel, const Momentum& ppr, const Limits& q2_range, const Limits& w_range);

    /// generate one event with unweighted photon and electron
    /** \note
       * 1) according to WWA:
       *    - transversal photonspectrum (\f$Q^2 \to 0\f$):
       *      \f$ P(y, Q^2) = \frac{\alpha}{2\pi}\frac{1}{Q^2 y}\left(2(1-y)\left(1-\frac{Q^2_{\rm min}}{Q^2}\right)+y^2\right),\f$
       *    - longitudinal photonspectrum (\f$Q^2 \to 0\f$):
       *      \f$ P(y, Q^2) = \frac{\alpha}{2\pi}\frac{1}{Q^2 y}2(1-y).\f$
       * 2) full transversal photonspectrum given by:
       *    - ABT, I. & J.R. SMITH (1992): MC upgrades to study untagged events. - H1-10/92-249.
       *    - SMITH, J.R. (1992): An experimentalist's guide to photon flux calculations. - H1-12/92-259.
       *    - SMITH, J.R. & B.D. BUROW (1994): Photon fluxes with beam mass effects and polarizations. - H1-01/94-338.
       * 3) full transversal and longitudinal spectrum by ABT&SMITH
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
    /// \return Helicity of the photon (+/-1: transverse photon, 0: longitudinal photon)
    short helicity(double longfr) const;
    /// alpha/2pi
    static const double ALPHARED;
    Mode mode_;
    Limits y_range_, q2_range_, w_range_, dy_range_;
    /// 5-vector of beam electron
    Momentum pel_;
    /// 5-vector of beam proton
    Momentum ppr_;
    double s_, w12_, elpr_, eel_;
    mutable std::array<unsigned int, 2> num_errors_;
    const std::array<unsigned int, 2> max_errors_;
    mutable double epa_max_;
  };
}  // namespace CepGen

#endif
