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
        if (_n_subtractions != _subtractions.size())
        {
            warning("form_factor", "Subtractions not set!");
            return NaN<double>();
        };

        for (int i = 0; i < _n_subtractions; i++)
        {
            subtraction_polynomial += _subtractions[i]*pow(s, i);
        };
        
        // If evaluating at the subtraction point just return the polynomial
        if (is_zero(s)) return subtraction_polynomial;

        return regular_piece(s);
    };

    // The discontinutiy which gets dispersed
    complex raw_form_factor::discontinuity(double s)
    {
        complex disc = 0;
        for (auto isobar : _direct_isobars)
        {
            disc += external_current(s)
                  * kinematic_factors(s)
                  * isobar->evaluate(s);
        };
        return disc;
    };
    
    // Non-singular part of the dispersion integral
    complex raw_form_factor::regular_piece(complex s)
    {
        using namespace boost::math::quadrature;

        // Only explciitly evaluate above cut, below is given by Schwarz
        if (imag(s) < 0) return conj(regular_piece(conj(s)));

        // Bounds of integration
        double low  = _kinematics->sth();
        double mid  = _settings._intermediate_energy;
        double high = _settings._cutoff;
        double eps  = _settings._infinitesimal;

        // If we're sufficiently far from the cut, evaluate normally
        bool far_from_cut = (real(s) <= _kinematics->sth()) || (imag(s) > eps);
        if  (far_from_cut)
        {
            auto fdx = [this,s](double x)
            {
                return discontinuity(x)/(x-s)*pow(s/x, _n_subtractions);
            };

            complex integral = gauss_kronrod<double,N_GAUSS_CAUCHY>::integrate(fdx, low, mid,  _settings._cauchy_integrator_depth, 1.E-9, NULL) 
                             + gauss_kronrod<double,N_GAUSS_CAUCHY>::integrate(fdx, mid, high, _settings._cauchy_integrator_depth, 1.E-9, NULL);
            return integral/PI;
        };

        // If we're close to the cut, use Cauchy trick with ieps
        double  rs     = real(s);
        complex disc_s = (rs <= _kinematics->sth()) ? 0 : discontinuity(rs);

        auto fdx = [this,rs,disc_s,eps,s](double x)
        {
            complex integrand;
            integrand  = discontinuity(x)-disc_s;
            integrand *= pow(rs/x, _n_subtractions);
            integrand /= (x-(rs+I*eps));
            return integrand;
        };

        // If using gauss gauss-legendre, split the integral into two pieces to avoid systematic errors at low energies
        complex integral  = gauss_kronrod<double,N_GAUSS_CAUCHY>::integrate(fdx, low, mid,  _settings._cauchy_integrator_depth, 1.E-9, NULL) 
                          + gauss_kronrod<double,N_GAUSS_CAUCHY>::integrate(fdx, mid, high, _settings._cauchy_integrator_depth, 1.E-9, NULL);
        complex logarithm = disc_s * log(1.-(s+I*eps)/low);
        
        return (integral-logarithm)/PI;
    };
}; // namespace iterateKT