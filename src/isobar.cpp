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
    complex raw_isobar::omnes(double s, complex ieps)
    {
        // We split the integration in two parts
        // to properly handle the Principle Value and ieps perscription

        double RHCs = (s > _kin->sth()) ? phase_shift(s) : 0.;
        auto fdx = [this, s, RHCs, ieps](double x)
        {
            complex integrand;
            integrand  = phase_shift(x) - RHCs;
            integrand *= (s/x); // One subtraction
            integrand /= (x - (s+ieps)); 
            return integrand;
        };

        // bounds of integration
        double low = _kin->sth();
        double high = std::numeric_limits<double>::infinity();

        complex logarithm, integral, result;
        integral  = boost::math::quadrature::gauss_kronrod<double,61>::integrate(fdx, low, high, _integrator_depth, 1.E-9, NULL);
        logarithm = are_equal(s, _kin->sth(), 1E-5) ? 0 : RHCs * log(1.-(s+ieps) / low);
        result = (integral-logarithm)/M_PI;
    
        return exp(result);
    };
};