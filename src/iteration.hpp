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
#include <TMath.h>
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
            _sth = kin->sth(); _pth = kin->pth(); _rth = kin->rth();

            // Load up the interpolators from the input data
            for (int i = 0; i < n; i++)
            {
                _re_inhom.push_back(new ROOT::Math::Interpolator(dat._s_list, dat._re_list[i]));
                _im_inhom.push_back(new ROOT::Math::Interpolator(dat._s_list, dat._im_list[i]));
            };

            // Also half-regularize them and save as an interpolations
            for (int i = 0; i < n; i++)
            {
                std::vector<double> re, im;
                for (auto s : dat._s_list)
                {
                    complex half_reg = half_regularized_integrand(i, s);
                    re.push_back( std::real(half_reg) );
                    im.push_back( std::imag(half_reg) );
                };
                _re_halfreg.push_back(new ROOT::Math::Interpolator(dat._s_list, re));
                _im_halfreg.push_back(new ROOT::Math::Interpolator(dat._s_list, im));
            };
            _initialized = true;
        };

        // Clean up manual pointers
        ~raw_iteration()
        {
            if (_zeroth) return;
            for (auto ptr : _re_inhom)   delete ptr;
            for (auto ptr : _im_inhom)   delete ptr;
            for (auto ptr : _re_halfreg) delete ptr;
            for (auto ptr : _im_halfreg) delete ptr;
        };

        // ----------------------------------------------------------------------- 
        // Evaluate the actual amplitude piece which usually involves a dispersion integral
        
        // This returns the function which multiplies the Omnes. Its given in the form
        // s^i + s^n * dispersion_integral(s)
        complex basis_factor(unsigned int i, complex x);

        // This is the inhomogeneity multiplied by kappa^N to remove all singularities
        complex ksf_inhomogeneity(unsigned int i, double x);

        // We take tilde_inhomogeneity and divide by nu^N which removes zeros at sth and rth
        // This is KSF and we only have the non-removable pth singularities to contend with
        complex half_regularized_integrand(unsigned int i, double x);
        complex regularized_integrand(unsigned int i, double x);

        // -----------------------------------------------------------------------
        private:

        // integration and interpolation settings
        settings _settings;

        // way to get thresholds and momenta
        kinematics _kinematics;
        
        // Whether this is the homogeneous solution with a trivial integral
        bool _zeroth = false;

        unsigned int _n_subtraction  = 1; // total subtractions (also # of basis functions)
        unsigned int _n_singularity  = 3; // degree of singular kinematic factors

        // These are convenient to save so to not have to keep redefining them
        double _sth, _pth, _rth;

        // The saved data and interpolation of the discontinuity
        std::vector<ROOT::Math::Interpolator*> _re_inhom, _im_inhom;

        // We're also going to save interpolations of the half-regularized disc
        std::vector<ROOT::Math::Interpolator*> _re_halfreg, _im_halfreg;

        // Flag to know if the interpolators are set up
        bool _initialized = false;

        // Calculate the expansion coefficients
        std::array<complex, 3> rthreshold_expansion(unsigned int i, double s, double eps);
        std::array<complex, 4> pthreshold_expansion(unsigned int i, double eps);

        // Analytic parts of the dispersion which we've subtracted away
        complex Q(int n, complex s, std::array<double,2> bounds);
        complex R(int n, complex s, std::array<double,2> bounds);

        // Consider two cases of dispersion integrals
        // They either contain the pth singularity or the cauchy one
        complex disperse_with_pth   (unsigned int i, complex s, std::array<double,2> bounds);
        complex disperse_with_cauchy(unsigned int i, complex s, std::array<double,2> bounds);
    };

}; // namespace iterateKT

#endif // ITERATION_HPP