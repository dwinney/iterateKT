// Top (abstract) level class which defines an X -> π transition form factor.
// We take an X -> 3π amplitude object and specify an isobar.
// This defines the type of current being calculated (i.e. P-wave isobar -> vector FF)
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Instituto de Ciencias Nucleares (ICN)
//               Universidad Autónoma Nacional de México (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#include "form_factor.hpp"

namespace iterateKT
{
    // Evaluate the dispersion relation, mimics strucutre of raw_iteration::integral
    complex raw_form_factor::evaluate(complex s)
    {
        // Subtraction polynomial
        complex subtraction_polynomial = 0;
        for (int i = 0; i < _n_subtractions; i++)
        {
            subtraction_polynomial += _subtractions[i]*pow(s, i);
        };
        
        // If evaluating at the subtraction point just return the polynomial
        if (is_zero(s)) return subtraction_polynomial;

        return 0.;
    };

    // Non-singular part of the dispersion integral
    complex raw_form_factor::regular_piece(complex s)
    {
        using namespace boost::math::quadrature;

        auto integrand = [this,s](double x)
        {
            complex disc = 0;
            for (auto isobar : _direct_isobars)
            {
                disc += external_current(x)
                      * kinematic_factors(x)
                      * isobar->evaluate(x);
            };
            return disc/(x-s)/pow(x, _n_subtractions);
        };

        return 0.;
         
        // // Integrate on either side of the pth singularity 
        // complex integral = gauss_kronrod<double,N_GAUSS_PSEUDO>::integrate(fdx, bounds[0], _pth, _settings._pseudo_integrator_depth, 1.E-9, NULL)
        //                  + gauss_kronrod<double,N_GAUSS_PSEUDO>::integrate(fdx, _pth, bounds[1], _settings._pseudo_integrator_depth, 1.E-9, NULL);
        // // Add back the analytic pieces we subtracted before
        // auto coeffs = (real(s) < _pth) ? _below_pth_expansion[i] : _above_pth_expansion[i];
        // for (int i = 0; i <= _l; i++) integral += coeffs[i] * Q(_n-2*i,s,bounds);
        // return integral;
    };
}; // namespace iterateKT