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

    inline iteration new_iteration(unsigned int nsub, std::array<std::vector<double>,3>& disc, settings sets)
    { 
        return std::make_shared<raw_iteration>(nsub, disc, sets); 
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
        // disc[0] = s
        // disc[1] = Re(discT)
        // disc[2] = Im(discT)
        raw_iteration(unsigned int n, std::array<std::vector<double>,3> disc, settings sets) 
        : _n(n), _zeroth(false), _settings(sets), _disc_data(disc)
        {
            _sth = disc[0].front(); _upper = disc[0].back();
            _re_disc.SetData(disc[0], disc[1]);
            _im_disc.SetData(disc[0], disc[2]);
        };

        // ----------------------------------------------------------------------- 
        // Evaluate the actual amplitude piece which usually involves a dispersion integral
        
        // This returns the function which multiplies the Omnes. Its given in the form
        // s^i + dispersion(s)
        complex basis_function(unsigned int i, complex x);

        // This is the discontinutiy which gets dispersed above
        complex discontinuity(double x);

        // -----------------------------------------------------------------------
        private:

        settings _settings;

        // threshold
        double _sth, _upper;

        // Number of basis functions
        unsigned int _n = 1;

        // Whether this is the homogeneous solution with a trivial integral
        bool _zeroth = false;

        // The saved data and interpolation of the discontinuity
        std::array<std::vector<double>,3> _disc_data;
        ROOT::Math::Interpolator _re_disc, _im_disc;
    };

}; // namespace iterateKT

#endif // ITERATION_HPP