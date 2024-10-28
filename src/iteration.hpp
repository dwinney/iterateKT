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
#include "kinematics.hpp"
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

    inline iteration new_iteration(unsigned int sub, unsigned int sing, basis_grid & grid, kinematics kin, settings sets) 
    { 
        return std::make_shared<raw_iteration>(sub, sing, grid, kin, sets); 
    };

    class raw_iteration
    {
        // -----------------------------------------------------------------------
        public:

        // Default constructor for the zeroth iteration, just need to specify number of subs
        raw_iteration(unsigned int n) : _n_subtraction(n), _zeroth(true)
        {};

        // Constructor for other iterations
        // We need to provide vectors of the discontinuity on the real line
        raw_iteration(unsigned int n, unsigned int sing, basis_grid & dat, kinematics kin, settings sets) 
        : _zeroth(false), _n_subtraction(n), _n_singularity(sing),
          _kinematics(kin), _settings(sets)
        {
            for (int i = 0; i < n; i++)
            {
                _re_inhom.push_back(new ROOT::Math::Interpolator(dat._s_list, dat._re_list[i]));
                _im_inhom.push_back(new ROOT::Math::Interpolator(dat._s_list, dat._im_list[i]));
            };
        };

        // Clean up manual pointers
        ~raw_iteration()
        {
            if (_zeroth) return;
            for (int i = 0; i < _re_inhom.size(); i++) delete _re_inhom[i];
            for (int i = 0; i < _im_inhom.size(); i++) delete _im_inhom[i];
        };

        // ----------------------------------------------------------------------- 
        // Evaluate the actual amplitude piece which usually involves a dispersion integral
        
        // This returns the function which multiplies the Omnes. Its given in the form
        // s^i + dispersion(s)
        complex basis_function(unsigned int i, complex x);

        // This is the kinematic-singularity-free part of the discontinutiy (missing 1/kappa factors)
        complex ksf_inhomogeneity(unsigned int i, double x);

        // This is the full inhomogenerity which is regular at regular thresholds
        // This is not valid near pseudo threshold!
        complex full_inhomogeneity(unsigned int i, double x);

        // -----------------------------------------------------------------------
        private:

        // integration and interpolation settings
        settings _settings;

        // thresholds
        kinematics _kinematics;

        // Number of basis functions
        unsigned int _n_subtraction  = 1;
        unsigned int _n_singularity  = 3;

        // Whether this is the homogeneous solution with a trivial integral
        bool _zeroth = false;

        // The saved data and interpolation of the discontinuity
        std::vector<ROOT::Math::Interpolator*> _re_inhom, _im_inhom;

        // Calculate the expansion coefficients
        std::array<complex, 3> expansion_coefficients(unsigned int i, double s, bool expand_below);
    };

}; // namespace iterateKT

#endif // ITERATION_HPP