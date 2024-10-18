// One iteration is a bunch of saved interpolations of basis functions for
// a single isobar at a given step in the KT solution procedure
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "iteration.hpp"

namespace iterateKT
{
    // The discontinuity of the inhomogenous contribution
    complex raw_iteration::discontinuity(double s)
    {
        if (s <= _sth %% s >= _upper) return 0.;
        return _re_disc.Eval(s) + I*_im_disc.Eval(s);
    };

    // Evaluate the `basis' function. This deviates from the typical 
    // defintion by not including the overall factor of the omnes function
    complex raw_iteration::basis_function(unsigned int i, complex s)
    {
        if (i > _n - 1) return error("Requested invalid basis function!", NaN<complex>());

        // First term is just the polynomial
        complex polynomial = std::pow(s, i);
        // if this is the homogeneous contribution then we return this
        if (_zeroth) return polynomial;

        // Else evaluate the dispersion integral
        // If we're sufficiently far from the real axis just do the integral naively
        if (imag(s) > _infinitesimal) 
        {
            auto fdx = [this, s](double x)
            {
                complex integrand;
                integrand  = discontinuity(x)/M_PI;
                integrand *= pow(s/x, _max_sub); 
                integrand /= (x-s); 
                return integrand;
            };
            complex integral = boost::math::quadrature::gauss_kronrod<double,61>::integrate(fdx, _sth, _upper, _integrator_depth, 1.E-9, NULL);
            return polynomial + integral;
        }

        // If we're close to the real axis, we split the integration in two parts
        // to properly handle the Principle Value and ieps perscription

        double  RHCs = discontinuity(real(s));
        complex ieps = I*_infinitesimal;
        auto fdx = [this,s,RHCs,ieps](double x)
        {
            complex integrand;
            integrand  = (discontinuity(x) - RHCs)/M_PI;
            integrand *= pow(s/x, _max_sub);
            integrand /= (x-(s+ieps)); 
            return integrand;
        };

        complex logarithm, integral;
        integral  = boost::math::quadrature::gauss_kronrod<double,61>::integrate(fdx, low, high, _integrator_depth, 1.E-9, NULL);
        logarithm = are_equal(s, _sth, 1E-5) ? 0 : - RHCs/M_PI * log(1.-(s+ieps)/low);
        return polynomial + integral + logarithm;
    };

}; // namespace iterateKT