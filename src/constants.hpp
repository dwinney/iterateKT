// Header file with global phyiscal constants.
// Everything is in GeV unless explicitly stated otherwise.
//
// ---------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <complex>
#include <limits>

namespace iterateKT
{
    using complex = std::complex<double>;

    // Number of gaussian quadrature points to use in nonadaptive integration
    const int N_GAUSS_DISPERSIVE = 200;
    const int N_GAUSS_OMNES      = 250;
    const int N_GAUSS_ANGULAR    = 100;

    // ---------------------------------------------------------------------------
    // Mathematical constants 

    #ifndef PI
        const double PI   = M_PI;
    #endif
    #ifndef I
        const complex I   (0., 1.);
    #endif

    const double DEG2RAD  = (M_PI / 180.);
    const double EPS      = 1.e-8;
    const complex IEPS    = I*EPS;

    // PDG Meson masses in GeV
    const double M_PION      = 0.13957000;
    const double M_KAON      = 0.49367700;
    const double M_ETA       = 0.54753;
    const double M_RHO       = 0.77526;
    const double M_OMEGA     = 0.78265;
    const double M_PHI       = 1.01956;
    const double M_F2        = 1.2754;
    const double M_B1        = 1.229;

    // ------------------------------------------------------------------------------
    // NaN's, 0, and 1 for throwing errors with custom data types

    template<typename T>
    inline T NaN()
    {
        return std::numeric_limits<T>::quiet_NaN();
    };

    template<>
    inline complex NaN()
    {
        return complex(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    };

    template<typename T>
    T zero();

    template<>
    inline complex zero() { return 0; };

    template<typename T> 
    T identity();

    template<>
    inline complex identity() { return 1; };

}; // iterateKT

#endif // CONSTANTS_HPP