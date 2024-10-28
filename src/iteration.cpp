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
    // -----------------------------------------------------------------------
    // The discontinuity of the inhomogenous contribution
    complex raw_iteration::ksf_inhomogeneity(unsigned int i, double s)
    {
        if (s <= _kinematics->sth() || s >= _settings._interp_energy_high) return 0.;
        return _re_inhom[i]->Eval(s) + I*_im_inhom[i]->Eval(s);
    };

    complex raw_iteration::full_inhomogeneity(unsigned int i, double s)
    {
        double Sth = _kinematics->sth();
        double Rth = _kinematics->rth();
        if (s < Sth) return NaN<double>();

        double   xi = _settings._matching_interval;
        complex nus = pow(_kinematics->nu(s,xi), _n_singularity);
        bool close_to_sth = are_equal(s, Sth, xi);
        bool close_to_rth = are_equal(s, Rth, xi);

        // if we're not near any threshold just divide like normal
        if (!(close_to_sth || close_to_rth)) return ksf_inhomogeneity(i, s)/nus;
        
        double eps = _settings._expansion_eps;
        double s_exp = (close_to_sth) ? Sth : Rth;
// 
        // The only time we expand below (i.e. s_exp - eps) is if s \simeq rth
        bool expand_below = (close_to_rth) && (s <= Rth);
        std::array<complex,3> ecs = expansion_coefficients(i, s_exp, expand_below);
        return (ecs[0]+ecs[1]*(s-s_exp)+ecs[2]*pow(s-s_exp,2))/nus;
    };

    std::array<complex,3> raw_iteration::expansion_coefficients(unsigned int i, double s, bool expand_below)
    {

        // These coefficients only depend on the order of the singularity
        int n = _n_singularity;
        std::array<int,4> as, bs, cs;
        switch (n)
        {
            case 1:  // S-waves
            {
                as = {+15, -12, +4, 8};
                bs = { -5,  +8, -4, 4};
                cs = { +3,  -4, +4, 8};
                break;
            };
            case 3:  // P-waves
            {
                as = {+35, -20, +4, 8};
                bs = {-21, +16, -4, 4};
                cs = {+15, -12, +4, 8};
                break;
            };
            case 5: // D-waves
            {
                as = {+63, -28, +4, 8};
                bs = {-45, +24, -4, 4};
                cs = {+35, -20, +4, 8};
                break;
            };
            default : 
           {
            warning("raw_iteration::expansion_coefficients", "Invalid singularity order (n="+std::to_string(n)+")!");
            return {NaN<complex>(), NaN<complex>(), NaN<complex>()};
           };
        };

        // Need up to second derivative
        double e = _settings._expansion_eps;
        if (expand_below) e *= -1;
        double s_exp = s+e;
        complex f    = _re_inhom[i]->Eval(s_exp)   + I*_im_inhom[i]->Eval(s_exp);
        complex fp   = _re_inhom[i]->Deriv(s_exp)  + I*_im_inhom[i]->Deriv(s_exp);
        complex fpp  = _re_inhom[i]->Deriv2(s_exp) + I*_im_inhom[i]->Deriv2(s_exp);

        complex a = (as[0]*f+as[1]*e*fp+as[2]*e*e*fpp)/(as[3]*pow(sqrt(std::abs(e)), n  ));
        complex b = (bs[0]*f+bs[1]*e*fp+bs[2]*e*e*fpp)/(bs[3]*pow(sqrt(std::abs(e)), n+2));
        complex c = (cs[0]*f+cs[1]*e*fp+cs[2]*e*e*fpp)/(cs[3]*pow(sqrt(std::abs(e)), n+4));

        if (expand_below) b *= -1;
        return {a,b,c};
    };

    // -----------------------------------------------------------------------
    // Evaluate the `basis' function. This deviates from the typical 
    // defintion by not including the overall factor of the omnes function
    complex raw_iteration::basis_function(unsigned int i, complex s)
    {
        if (i > _n_subtraction - 1) return error("Requested invalid basis function (i = "+std::to_string(i)+", n = " +std::to_string(_n_singularity)+")!", NaN<complex>());

        // First term is just the polynomial
        complex polynomial = std::pow(s, i);
        // if this is the homogeneous contribution then we return this
        if (_zeroth) return polynomial;

        return NaN<complex>();
        // // Else evaluate the dispersion integral
        // // If we're sufficiently far from the real axis just do the integral naively
        // double eps = _settings._infinitesimal;
        // if (imag(s) > eps) 
        // {
        //     auto fdx = [this,i,s](double x)
        //     {
        //         complex integrand;
        //         integrand  = ksf_inhomogeneity(i, x)/M_PI;
        //         // integrand *= pow(s/x, _n); 
        //         integrand /= (x-s); 
        //         return integrand;
        //     };
        //     complex integral = boost::math::quadrature::gauss_kronrod<double,61>::integrate(fdx, _sth, _upper, _settings._dispersion_integrator_depth, 1.E-9, NULL);
        //     return polynomial + integral;
        // }

        // // If we're close to the real axis, we split the integration in two parts
        // // to properly handle the Principle Value and ieps perscription

        // complex RHCs = ksf_inhomogeneity(i, real(s));
        // auto fdx = [this,s,RHCs,eps,i](double x)
        // {
        //     complex integrand;
        //     integrand  = (ksf_inhomogeneity(i, x) - RHCs)/M_PI;
        //     // integrand *= pow(s/x, _n);
        //     integrand /= (x-(s+I*eps)); 
        //     return integrand;
        // };

        // complex logarithm, integral;
        // integral  = boost::math::quadrature::gauss_kronrod<double,61>::integrate(fdx, _sth, _upper, _settings._dispersion_integrator_depth, 1.E-9, NULL);
        // logarithm = are_equal(s, _sth, 1E-5) ? 0 : - RHCs/M_PI * log(1.-(s+I*eps)/_sth);
        // return polynomial + integral + logarithm;
    };

}; // namespace iterateKT