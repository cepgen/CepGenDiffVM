#ifndef CepGenDiffVM_EPA_h
#define CepGenDiffVM_EPA_h

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/Constants.h"
#include <iostream>

namespace CepGen
{
  /// \author Benno List
  /// \author Thomas Jansen
  /// \date 1993-1995
  class EPA
  {
    public:
      /// Photon generation mode
      enum struct Mode
      {
        wwa = 1, ///< WWA approximation (including e-mass effect and longitudinal flux). Recommended
        transverse = 2, ///< transverse spectrum
        transverse_longitudinal = 3, ///< transverse & longitudinal spectrum
        transverse_longitudinal_pframe = 4 ///< as 3, but flux in p rest frame
      };
      friend std::ostream& operator<<( std::ostream& os, const Mode& mode ) {
        switch ( mode ) {
          case Mode::wwa: return os << "WWA";
          case Mode::transverse: return os << "transverse only";
          case Mode::transverse_longitudinal: return os << "transverse+longitud.";
          case Mode::transverse_longitudinal_pframe: return os << "transverse+longitud.[p frame]";
        }
        return os;
      }
      EPA( const Mode& mode = Mode::wwa ) :
        mode_( mode ), y_range_( { 0., 1. } ), dy_range_( { 0., 1. } ),
        s_( 0. ), w12_( 0. ), elpr_( 0. ), eel_( 0. ),
        num_errors_( { { 0, 0, 0 } } ), epa_max_( -1. ) {}
      struct Result
      {
        Result() : valid( false ), q2( 0. ), heli( 0 ) {}
        bool valid;
        /// 5-vector of photon in H1 lab. system (5th component is sign (q2)*sqrt (abs (q2)))
        Particle::Momentum pph;
        /// 5-vector of scattered electron
        Particle::Momentum ppe;
        /// Virtuality of photon (positive!): Q2 = -q2
        double q2;
        /// Photon helicity:
        /// - 0: longitudinal,
        /// - +/-1: transverse polarization
        short heli;
      };

      /// Initialize histograms, constants, kinematical bounds
      void init( const Particle::Momentum& pel, const Particle::Momentum& ppr, const Limits& q2_range, const Limits& w_range ) {
        pel_ = pel;
        ppr_ = ppr;
        q2_range_ = q2_range;
        w_range_ = w_range;

        CG_INFO( "EPA:init" ) << "EPA module initialised.\n\t"
          << "mode: " << mode_ << "\n\t"
          << "beams momenta:\n\t  " << pel_ << "\n\t  " << ppr_ << "\n\t"
          << "W range: " << w_range_ << "\n\t"
          << "Q² range: " << q2_range_;

        //--- calculate CMS s=(P+K)²
        elpr_ = ppr_.energy()*pel_.energy()+ppr_.p()*pel_.p();
        s_ = 2.*elpr_+pel_.mass2()+ppr_.mass2();
        if ( mode_ == Mode::transverse_longitudinal_pframe )
          //--- evaluate photon flux in proton rest frame:
          // set EEL to approx. 50TeV
          eel_ = elpr_/ppr_.mass();
        else
          eel_ = pel_.energy();
        CG_DEBUG( "EPA:init" ) << "S = " << s_ << ", EEL = " << eel_ << ".";
        const double w_min = w_range_.min(), w_min2 = w_min*w_min,
                     w_max = w_range_.max();
        const double q2_min = q2_range_.min(), q2_max = q2_range_.max();
        w12_ = w_min2-ppr_.mass2();

        //--- calculate Y - bounds from
        // ALI, A. et al. (1987): Heavy quark physics at HERA. - Proc. HERA workshop, Hamburg 1987 (ed. R.D. PECCEI), 395-494.
        const double y_sqr = sqrt( pow( s_-w12_, 2 )-4.*w12_*pel_.mass2() );
        dy_range_.max() = 0.5*( s_+w12_+y_sqr )/( s_+pel_.mass2() );
        //---  use trick for quadratic equations; see
        // W.H. PRESS et al. (1988): Numerical Recipes in C. Cambridge (Cambridge Univ. Press), p. 156.
        dy_range_.min() = std::max( w12_/( dy_range_.max()*( s_+pel_.mass2() ) ), y_range_.min() );
        //--- calculate absolute maximum of y, irrespective of final state
        dy_range_.max() = std::min( std::min( s_/( s_+pel_.mass2() ), y_range_.max() ), 0.5*( w_max*w_max-ppr_.mass2()+q2_max )/elpr_ );
        CG_INFO( "EPA:init" )
          << "S / (S + DME**2) = " << ( s_/( s_+pel_.mass2() ) ) << "\n\t"
          << "y range: " << y_range_ << "\n\t"
          << "(W(max)²-mp²+Q²(max))/(2*ELPR) = " << ( 0.5*( w_max*w_max-ppr_.mass2()+q2_max )/elpr_ ) << "\n\t"
          << "dy range: " << dy_range_ << ".";
        const double dy_min = dy_range_.min(), dy_max = dy_range_.max();

        //--- set max. photon weight for efficient rejection plane
        const double qg2_min = std::max( pel_.mass2()*pow( dy_min, 2 )/( 1.-dy_min ), q2_min );
        const double qg2_max = std::min( dy_max*s_, q2_max );
        Limits qg2_range( qg2_min, qg2_max );
        CG_DEBUG( "EPA:init" )
          << "QG2 range: " <<  qg2_range << ".";

        if ( mode_ == Mode::wwa ) //--- WWA - approximation
          epa_max_ = ALPHARED * ( 4.*( 1.-dy_min ) + pow( dy_min, 2 ) );
        else { //--- full transversal spectrum (2) or full longitudinal and transversal (3) spectrum
          const double eqe = qg2_min/eel_/eel_;
          const double emqe2 = pow( dy_min-0.25*eqe, 2 );
          const double emsqr = ( pow( dy_min*elpr_, 2 )+qg2_min*ppr_.mass2() )/( elpr_*elpr_+pel_.mass2()*ppr_.mass2() );

          if ( emsqr < 0. )
            CG_FATAL( "EPA:init" ) << "problem with sqrt(emsqr): "
              << emsqr << " at EPAMAX determination.";

          epa_max_ = ALPHARED * dy_min * sqrt( emsqr ) / ( emqe2+eqe );
          if ( mode_ == Mode::transverse )
            epa_max_ *= ( 2.*( 1.-dy_min )+emqe2+eqe );
          else
            epa_max_ *= ( 4.*( 1.-dy_min )+emqe2+eqe );
        }
        epa_max_ *= log( dy_max/dy_min ) * log( qg2_max/qg2_min );

        CG_DEBUG( "EPA:init" ) << "maximal EPA: " << epa_max_;
      }

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
       *
       *       kinematical bounds \f$Y_{\rm min}\f$,  \f$Y_{\rm max} (W_{\rm min})\f$, \f$Q^2_{\rm min}\f$, \f$Q^2_{\rm max} (Q^2_{\rm cutoff})\f$.
       */
      Result operator()( double x1, double x2 ) const {
        // GEPHOT (PEL, PPR, PPH, PPE, Q2, HELI)
        Result out;

        //--------------------------------------------------------------
        // begin main loop over Y,Q2 rnd prod.
        //--------------------------------------------------------------

        const double dy_min = dy_range_.min(), dy_max = dy_range_.max();

        //--- produce Y spect. ( 1/y weighted shape )
        const double y = dy_min*pow( dy_max/dy_min, x1 );
        //--- calculate actual Q2_min, Q2_max from Y
        const double gq2_min = std::max( pel_.mass2()*y*y/( 1.-y ), q2_range_.min() );
        const double gq2_max = std::min( y*s_, q2_range_.max() );
        const Limits gq2_range( gq2_min, gq2_max );
        //--- take Q2_cut from steering, if it is kinematicly reachable.
        if ( !gq2_range.valid() )
          return out;

        CG_DEBUG( "EPA:get" )
          << "Y, Q2min, Q2max: " << y << ", " << gq2_range << ".";

        double epa = 0., epa_t = 0., epa_l = 0., lf = 0.;

        //--- produce Q2 spect. (1/x weighted shape )
        const double q2 = gq2_min*pow( gq2_max/gq2_min, x2 );

        //--------------------------------------------------------------
        // EPA - WWA spectrum
        //--------------------------------------------------------------

        //--- calc. photon weight
        if ( mode_ == Mode::wwa ) {
          //--- WWA - approximation
           const double r = ALPHARED/( y*q2 );
           epa_t = r*( 2.*( 1.-y )*( 1.-pel_.mass2()*y*y/( ( 1.-y )*q2 ) )+y*y );
           epa_t = r*( 2.*( 1.-y ) );
           epa = epa_t+epa_l;
           lf = epa_l/epa;
        }
        else {
          //------------------------------------------------------------
          // full transversal spectrum (2) or full longitudinal and
          // transversal (3) spectrum from:
          // * ABT, I. & J.R. SMITH (1992): MC Upgrades to Study
          //    Untagged Events. - H1-10/92-249.
          // See also:
          // * SMITH, J.R. (1992): An Experimentalist's Guide to Photon
          //    Flux Calculations. - H1-12/92-259.
          // * SMITH, J.R. (1993): Polarization Decomposition of Fluxes
          //    and Kinematics in ep Reactions. - H1-04/93-282.
          //------------------------------------------------------------
          const double eqe = q2/eel_/eel_;
          const double emqe2 = pow( y-0.25*eqe, 2 );
          const double emsqr = ( pow( y*elpr_, 2 )+q2*ppr_.mass2() )/( elpr_*elpr_+pel_.mass2()*ppr_.mass2() );

          if ( emsqr < 0. ) {
            CG_WARNING( "EPA:get" )
              << "problem with sqrt(emsqr) = " << emsqr << ": "
              << "y, Q2 pair rejected";
            if ( ++num_errors_[0] > 10 )
              throw CG_FATAL( "EPA:get" ) << "too many sqrt problems: try WWA!";
          }

          const double r = ALPHARED/q2*sqrt( emsqr )/( emqe2+eqe );
          epa_t = r*( 2.*( 1.-y )+emqe2+eqe );
          epa_l = 0.;
          if ( mode_ != Mode::transverse ) //--- longitudinal & transversal spectrum
            epa_l = r*( 2.*( 1.-y ) );
        }

        epa = epa_t+epa_l;
        lf = epa_l/epa;

        //--- unweight MC
        double r = y*q2*log( dy_max/dy_min )*log( gq2_max/gq2_min );
        const double w = sqrt( y*2.*elpr_-q2+ppr_.mass2() );
        //--- check if W_min < W < W_max, else reject photon
        if ( !w_range_.passes( w ) )
          r = 0.;
        epa *= r;
        epa_t *= r;
        epa_l *= r;

        //--- update upper EPA bound
        if ( epa > epa_max_ ) {
          if ( epa > 1.1*epa_max_ )
            CG_WARNING( "EPA:get" ) << "S: EPA > 1.1*EPAMAX!";
          else if ( epa > 1.01*epa_max_ )
            CG_WARNING( "EPA:get" ) << "W: EPA > 1.01*EPAMAX!";
          else
            CG_WARNING( "EPA:get" ) << "I: EPA > EPAMAX!";
          epa_max_ = epa;
          CG_INFO( "EPA:get" ) << "update of maximal weight: " << epa_max_ << ".";
        }

        CG_DEBUG_LOOP( "EPA:get" )
          << "Y: " << y << ", Q²: " << q2 << "\n\t"
          << "EPA(T): " << epa_t << ", EPA(L): " << epa_l << "\n\t"
          << "EPA: " << epa << ", long.fraction: " << lf << "\n\t"
          << "GQ2 range: " << gq2_range << ".";

        //----> end rnd loop: rejection method
        if ( rand()*1./RAND_MAX*epa_max_ > epa )
          return out;
        //--- continue with unweighted kin.

        //--- scattering angle of electron in LAB.:
        //  E.Lohrmann DESY HERA 83/08
        // x = Q2 / (y s)
        // E_sc = E_e(1-y) + E_p x y
        // cos t = [E_e(1-y) - E_p x y] / E_sc
        const double emy = pel_.energy()*( 1.-y ), exy = ppr_.energy()*q2/s_;
        const double eesc  = emy+exy;
        const double cthe = ( emy-exy )/eesc, sthe = 2.*sqrt( emy-exy )/eesc;

        //--- control scattering angle
        if ( fabs( cthe ) > 1. || fabs( sthe ) > 1. ) {
          CG_DEBUG_LOOP( "EPA:get" ) << "theta of electron: "
            << asin( sthe ) << "/" << acos( cthe ) << ", "
            << "reject event for Y: " << y << ", Q2: " << q2 << ".";
          if ( ++num_errors_[1] > 500 )
            throw CG_FATAL( "EPA:get" ) << "too many problems for (S/C)THE!";
          //--- new kinematics
          return out;
        }

        out.valid = true;
        const double phi = 2.*M_PI*rand()/RAND_MAX;

        //--- 5-vector of electron in LAB system
        const double pesc = -sqrt( eesc*eesc-pel_.mass2() );
        out.ppe = Particle::Momentum::fromPThetaPhi( pesc, acos( cthe ), phi, eesc );
        out.ppe.setMass( pel_.mass() );

        //--- 5-vector of photon k - k'
        out.pph = pel_-out.ppe;
        out.pph.setMass( -sqrt( q2 ) );

        //--- determine helicity of the photon
        out.heli = helicity( lf );
        return out;
      }

    private:
      /// Generate the helicity of a photon
      /// \author Benno List
      /// \date 27 May 1993
      /// \param[in] longfr Fraction of longitudinally polarized photons
      /// \note IHELI in DIFFVM
      /// \return Helicity of the photon (+/-1: transverse photon, 0: longitudinal photon)
      short helicity( double longfr ) const {
        if ( rand()/RAND_MAX < longfr )
          return 0;
        if ( rand()/0.5 )
          return +1;
        return -1;
      }
      /// alpha/2pi
      static constexpr double ALPHARED = Constants::alphaEM*0.5*M_1_PI;
      Mode mode_;
      Limits y_range_, q2_range_, w_range_, dy_range_;
      /// 5-vector of beam electron
      Particle::Momentum pel_;
      /// 5-vector of beam proton
      Particle::Momentum ppr_;
      double s_, w12_, elpr_, eel_;
      mutable std::array<unsigned int,3> num_errors_;
      mutable double epa_max_;
  };
}

#endif
