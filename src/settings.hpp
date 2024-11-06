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

        // We have three types of integrations, we can choose whether to evaluate them
        // using adaptive integration or not
        bool _adaptive_omnes       = false; 
        bool _adaptive_dispersion  = false; 
        bool _adaptive_angular     = false; 

        // Number of subdivisions for adaptive integrator 
        // These are only looked at if the appropriate flag above is true
        int  _omnes_depth        = 15;    
        int  _dispersion_depth   = 15;
        int  _angular_depth      = 15;   

        // The infinitesimal to use for ieps inside Cauchy kernels
        double _infinitesimal      = 1E-5;

        // Interval +- regular thresholds around which to remove singularities
        std::array<double,3> _matching_intervals = {0.05, 0.05, 0.05};
        std::array<double,3> _expansion_offsets  = {0.05, 0.05, 0.05};

        // Interpolation settings
        std::array<int,3> _interpolation_points = {100, 40, 100};

        double _intermediate_energy  = 5;    // interpolate from sth to this value
        double _cutoff               = 1000; // then from _interp_energy_low to _interp_energy_high 
        double _interpolation_offset = 1;    // Amount to offset the middle point between interpolation regions        
    };
}; // namespace iterateKT

#endif // SETTINGS_HPP