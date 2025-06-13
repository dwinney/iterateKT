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
#include <Math/Interpolator.h>

namespace iterateKT
{
    struct settings
    {
        settings(){};

        // We use variable iterations to improve convergence by articially surpressing KT effects
        // At each next iteration we decrease the surpression linearly until we arrive back 
        // at the original KT equations
        double _iteration_rate_intercept  = 1.0;
        double _iteration_rate_slope      = 0.;

        // Number of subdivisions for adaptive integrator 
        // These are only looked at if the appropriate flag above is true
        double _omnes_cutoff     = std::numeric_limits<double>::infinity();
        int  _omnes_integrator_depth        = 0;    
        int  _cauchy_integrator_depth       = 0;
        int  _pseudo_integrator_depth       = 0;
        int  _angular_integrator_depth      = 0;   

        // Which type of interpolation to use, kCSPLINE or kAKIMA
         ROOT::Math::Interpolation::Type _interpolation_type = ROOT::Math::Interpolation::Type::kCSPLINE;

        // The infinitesimal to use for ieps inside Cauchy kernels
        double _infinitesimal      = 1E-5;
        double _derivative_h       = 1E-3;

        // Interval +- regular thresholds around which to remove singularities
        std::array<double,3> _matching_intervals = {0.05, 0.05, 0.05};
        std::array<double,3> _expansion_offsets  = {0.05, 0.05, 0.05};

        // We exclude a specific range around pth when evaluating the basis functions
        // and allow an interpolation to fit in the middle.
        int _exclusion_points = 10;
        std::array<double,2> _exclusion_offsets  = {0, 0};

        // Interpolation settings
        double _intermediate_energy  = 5;    // interpolate from sth to this value
        double _cutoff               = 1000; // then from _interp_energy_low to _interp_energy_high 
        double _interpolation_offset = 1;    // Amount to offset the middle point between interpolation regions        
        std::array<int,3> _interpolation_points = {100, 40, 100}; // Number of points to interpolate in each of the above regions
    };
}; // namespace iterateKT

#endif // SETTINGS_HPP