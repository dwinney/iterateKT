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
    constexpr int N_GAUSS_CAUCHY  = 151;
    constexpr int N_GAUSS_PSEUDO  = 61;
    constexpr int N_GAUSS_OMNES   = 61;
    constexpr int N_GAUSS_ANGULAR = 31;

    // ---------------------------------------------------------------------------
    // PDG Meson masses in GeV

    const double M_PION      = 0.13957000;
    const double M_KAON      = 0.49367700;

    // Different masses for isospin projectionss
    const double M_KAON_PM  = 0.493677;
    const double M_KAON_0   = 0.497611;
    const double M_KAON_AVG = (M_KAON_PM + M_KAON_0)/2;
    const double M_PION_PM  = 0.13957039;
    const double M_PION_0   = 0.1349768;
    const double M_PION_AVG = (M_PION_PM + M_PION_0)/2;

    const double M_ETA       = 0.54753;
    const double M_RHO       = 0.77526;
    const double M_OMEGA     = 0.78265;
    const double M_PHI       = 1.01956;
    const double M_F2        = 1.2754;
    const double M_B1        = 1.229;
    const double M_JPSI      = 3.096;
}; // iterateKT

#endif // CONSTANTS_HPP