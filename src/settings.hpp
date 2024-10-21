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

        int    _integrator_depth   = 15;
        double _infinitesimal      = 1E-5;

        // Interpolation settings
        double _interp_energy_low  = 5;    // interpolate from sth to this value
        int    _interp_points_low  = 200;  // interpolate the above interval with this many points
        double _interp_energy_high = 1000; // then from _interp_energy_low to _interp_energy_high 
        int    _interp_points_high = 200;  // with this many points
    };
}; // namespace iterateKT

#endif // SETTINGS_HPP