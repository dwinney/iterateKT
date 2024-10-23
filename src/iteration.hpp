// One iteration is a bunch of saved interpolations of basis functions for
// a single isobar at a given step in the KT solution procedure
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef ITERATION_HPP
#define ITERATION_HPP

#include <memory>
#include <Math/Interpolator.h>
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "basis_grid.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace iterateKT
{
    // Forward declare for the typedef below
    class raw_iteration;

    // Define isobars only as pointers
    using iteration  = std::shared_ptr<raw_iteration>;

    // This function serves as our "constructor"
    // iterations only need to know the number of basis functions to consider
    inline iteration new_iteration(unsigned int nsub)
    { 
        return std::make_shared<raw_iteration>(nsub); 
    };

    inline iteration new_iteration(unsigned int n, basis_grid & grid, settings sets) 
    { 
        return std::make_shared<raw_iteration>(n, grid, sets); 
    };

    class raw_iteration
    {
        // -----------------------------------------------------------------------
        public:

        // Default constructor for the zeroth iteration, just need to specify number of subs
        raw_iteration(unsigned int n) : _zeroth(true)
        {};

        // Constructor for other iterations
        // We need to provide vectors of the discontinuity on the real line
        raw_iteration(unsigned int n, basis_grid dat, settings sets) 
        : _n(n), _zeroth(false), _settings(sets)
        {
            _sth = dat._s_list.front(); _upper = dat._s_list.back();

            for (int i = 0; i < n; i++)
            {
                _re_disc.push_back(new ROOT::Math::Interpolator(dat._s_list, dat._re_list[i]));
                _im_disc.push_back(new ROOT::Math::Interpolator(dat._s_list, dat._im_list[i]));
            };
        };

        // Clean up manual pointers
        ~raw_iteration()
        {
            if (_zeroth) return;
            for (int i = 0; i < _n; i++)
            {
                delete _re_disc[i];
                delete _im_disc[i];
            };
        };

        // ----------------------------------------------------------------------- 
        // Evaluate the actual amplitude piece which usually involves a dispersion integral
        
        // This returns the function which multiplies the Omnes. Its given in the form
        // s^i + dispersion(s)
        complex basis_function(unsigned int i, complex x);

        // This is the kinematic-singularity-free part of the discontinutiy which gets dispersed 
        // over. 
        complex ksf_discontinuity(unsigned int i, double x);

        // -----------------------------------------------------------------------
        private:

        // integration and interpolation settings
        settings _settings;

        // thresholds
        double _sth, _upper;

        // Number of basis functions
        unsigned int _n = 1;

        // Whether this is the homogeneous solution with a trivial integral
        bool _zeroth = false;

        // The saved data and interpolation of the discontinuity
        std::vector<ROOT::Math::Interpolator*> _re_disc, _im_disc;
    };

}; // namespace iterateKT

#endif // ITERATION_HPP