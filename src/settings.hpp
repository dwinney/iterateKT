// Simple container class with all the different settings for isobar evaluation
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include "utilities.hpp"

namespace iterateKT
{
    struct settings
    {
        settings(){};

        // Number of subdivisions for adaptive integrator 
        // For dispersion integrals and the pinocchio integral
        int    _dispersion_integrator_depth = 15;
        int    _angular_integrator_depth    = 15;

        // The infinitesimal to use for ieps inside Cauchy kernels
        double _infinitesimal      = 1E-5;

        // Interval +- regular thresholds around which to remove singularities
        std::array<double,3> _matching_intervals = {0.05, 0.05, 0.05};
        std::array<double,3> _expansion_offsets  = {0.05, 0.05, 0.05};

        // Interpolation settings
        double _interp_energy_low  = 5;    // interpolate from sth to this value
        int    _interp_points_low  = 100;  // interpolate the above interval with this many points
        double _interp_energy_high = 1000; // then from _interp_energy_low to _interp_energy_high 
        int    _interp_points_high = 100;  // with this many points

        // Amount to offset the middle point between interpolation regions
        double _interp_offset = 1;
    };
}; // namespace iterateKT

#endif // SETTINGS_HPP