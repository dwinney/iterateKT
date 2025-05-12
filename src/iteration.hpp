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
#include "basis.hpp"
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
    inline iteration new_iteration(){ return std::make_shared<raw_iteration>(); };

    inline iteration new_iteration(kinematics kin, basis_grid & grid, settings sets) 
    { 
        return std::make_shared<raw_iteration>(kin, grid, sets); 
    };

    class raw_iteration
    {
        // -----------------------------------------------------------------------
        public:

        // Default constructor for the zeroth iteration, just need to specify number of subs
        raw_iteration() : _zeroth(true)
        {};

        // Constructor for other iterations
        // We need to provide vectors of the discontinuity on the real line
        raw_iteration(kinematics kin, basis_grid & dat, settings sets) 
        : _zeroth(false),   _n(dat._n_singularity),
          _kinematics(kin), _settings(sets)
        { 
            _l = (_n-1)/2;
            initialize(dat);
        };

        // Clean up manual pointers
        ~raw_iteration()
        {
            if (_zeroth) return;
            for (auto ptr : _re)                delete ptr;
            for (auto ptr : _im)                delete ptr;
            for (auto ptr : _re_excluded_above) delete ptr;
            for (auto ptr : _im_excluded_above) delete ptr;
            for (auto ptr : _re_excluded_below) delete ptr;
            for (auto ptr : _im_excluded_below) delete ptr;
        };

        // ----------------------------------------------------------------------- 
        // Evaluate the actual amplitude piece which usually involves a dispersion integral
        
        // This returns the inhomogeneity integral which multiplies the omnes function
        complex integral(uint i, complex s);

        // This is the inhomogeneity multiplied by kappa^N to remove all singularities
        complex ksf_inhomogeneity(uint i, double x);

        // We take tilde_inhomogeneity and divide by nu^N which removes zeros at sth and rth
        // This is KSF and we only have the non-removable pth singularities to contend with
        complex half_regularized_integrand(uint i, double x);
        complex regularized_integrand(uint i, double x);

        // -----------------------------------------------------------------------
        private:

        // integration and interpolation settings
        settings _settings;

        // way to get thresholds and momenta
        kinematics _kinematics;

        // Info regarding subtractions
        subtractions _subtractions;

        // Save all the interpolations and everything
        void initialize(basis_grid & data);
        
        // Whether this is the homogeneous solution with a trivial integral
        bool _zeroth = false;

        uint _l  = 1; // Angular momentum 
        uint _n  = 3; // degree of singular kinematic factors

        // These are convenient to save so to not have to keep redefining them
        double _sth, _pth, _rth;

        // The threshold square root continuation
        inline complex k(complex s)
        { 
            if (abs(imag(s)) >= _settings._infinitesimal) return csqrt(_pth - s);
            return  (real(s) <= _pth) ? csqrt(_pth-s) : +I*csqrt(s-_pth);  
        };

        // The saved data and interpolation of the discontinuity
        std::vector<ROOT::Math::Interpolator*> _re, _im;

        // Finally, for each basis function we have a small region around pth which we interpolate
        std::vector<ROOT::Math::Interpolator*> _re_excluded_above, _im_excluded_above, 
                                               _re_excluded_below, _im_excluded_below;

        // Flag to know if the interpolators are set up
        bool _initialized = false;
        std::vector<std::array<complex,4>> _below_pth_expansion, _above_pth_expansion;
        std::vector<std::array<complex,3>> _sth_expansion, _below_rth_expansion, _above_rth_expansion;
        // Calculate the expansion coefficients
        std::array<complex, 3> rthreshold_expansion(uint i, double s, double eps);
        std::array<complex, 4> pthreshold_expansion(uint i, double eps);

        // Analytic parts of the dispersion which we've subtracted away
        complex Q(int n, complex s, std::array<double,2> bounds);

        // Consider two cases of dispersion integrals
        // They either contain the pth singularity or the cauchy one
        complex disperse_with_pth   (uint i, complex s, std::array<double,2> bounds);
        complex disperse_with_cauchy(uint i, complex s, std::array<double,2> bounds);
    };

}; // namespace iterateKT

#endif // ITERATION_HPP