// The decay amplitude is decomposed into terms of one-variable functions
// These are given by the isobar class below
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "isobar.hpp"

namespace iterateKT
{
    // ----------------------------------------------------------------------- 
    // Evaluate the Omnes functions from the given phase_shift
    complex raw_isobar::omnes(complex s)
    {
        // Only explciitly evaluate above cut, below is given by Schwarz
        if (imag(s) < 0) return std::conj(omnes(std::conj(s)));

        // bounds of integration
        double low = _kinematics->sth();
        double high = std::numeric_limits<double>::infinity();
        
        // If we're sufficiently far from the real axis just do the integral naively
        double eps = _settings._infinitesimal;
        if (imag(s) > eps) 
        {
            auto fdx = [this, s](double x)
            {
                complex integrand;
                integrand  = phase_shift(x);
                integrand *= (s/x); // One subtraction
                integrand /= (x-s); 
                return integrand;
            };
            complex integral = boost::math::quadrature::gauss_kronrod<double,61>::integrate(fdx, low, high, _settings._integrator_depth, 1.E-9, NULL);
            return exp(integral/M_PI);
        }

        // If we're close to the real axis, we split the integration in two parts
        // to properly handle the Principle Value and ieps perscription

        double  RHCs = phase_shift(real(s));
        auto fdx = [this,s,RHCs,eps](double x)
        {
            complex integrand;
            integrand  = phase_shift(x) - RHCs;
            integrand *= (s/x); // One subtraction
            integrand /= (x-(s+I*eps)); 
            return integrand;
        };

        complex logarithm, integral, result;
        integral  = boost::math::quadrature::gauss_kronrod<double,61>::integrate(fdx, low, high, _settings._integrator_depth, 1.E-9, NULL);
        logarithm = are_equal(s, _kinematics->sth(), 1E-5) ? 0 : RHCs * log(1.-(s+I*eps) / low);
        result = (integral-logarithm)/M_PI;
    
        return exp(result);
    };

    // ----------------------------------------------------------------------- 
    // Specify a given iteration to use when outputting the basis_function
    complex raw_isobar::basis_function(unsigned int iter_id, unsigned int basis_id, complex x)
    {
        if (iter_id  > _iterations.size()) return error("Requested iteration does not exist!", NaN<complex>());
        if (basis_id > _max_sub - 1)  return error("Requested basis function does not exist!", NaN<complex>());

        return omnes(x)*_iterations[iter_id]->basis_function(basis_id, x);
    };

    // Without an iter_id we just take the latest iteration
    complex raw_isobar::basis_function(unsigned int basis_id, complex x)
    { 
        return basis_function(_iterations.size()-1, basis_id, x); 
    };

    // ----------------------------------------------------------------------- 
    // Combine all the basis_functions to evaluate the isobar
    complex raw_isobar::evaluate(unsigned int iter_id, complex s)
    {
        if (iter_id  > _iterations.size()) return error("Requested iteration does not exist!", NaN<complex>());
        
        complex result = 0.;
        for (int i = 0; i < _max_sub; i++) result += _subtraction_coeffs[i]*basis_function(iter_id, i, s);
        return result;
    };
    complex raw_isobar::evaluate(complex s)
    {
        return evaluate(_iterations.size()-1, s);
    };


    // ----------------------------------------------------------------------- 
    // Take the saved interpolation settings and output the necessary arrays
    
    std::array<std::vector<double>,3> raw_isobar::calculate_next(std::vector<isobar> previous)
    {
        // Interpolation in three different regions
        double low  = _kinematics->sth();
        double mid  = _settings._interp_energy_low;
        double high = _settings._interp_energy_high;

        int N_low  = _settings._interp_points_low;
        int N_high = _settings._interp_points_high;

        // First region is [low, mid]
        for (int i = 0; i < N_low; i++)
        {
            double s = low + (mid - high)*double(i)/double(N_low-1);
        };

        return {};
    };

};