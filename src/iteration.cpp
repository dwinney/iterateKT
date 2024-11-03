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
    // We were fed this from the constructor so we just access the interpolation
    complex raw_iteration::ksf_inhomogeneity(unsigned int i, double s)
    {
        if (s <= _sth || s >= _cutoff) return 0.;
        return _re_inhom[i]->Eval(s) + I*_im_inhom[i]->Eval(s);
    };

    // -----------------------------------------------------------------------
    // Take the above inhomogeneity and divide by the appropriate number of 
    // nu given by n_singularity
    complex raw_iteration::half_regularized_integrand(unsigned int i, double s)
    {
        if (s < _sth)     return NaN<double>();
        if (_initialized) return _re_halfreg[i]->Eval(s) + I*_im_halfreg[i]->Eval(s);

        // If no interpolation is saved yet, calculate divide by nu explicitly
        double   xi = _settings._matching_interval;
        complex nus = pow(_kinematics->nu(s,xi), _n_singularity);
        bool close_to_sth = are_equal(s, _sth, xi);
        bool close_to_rth = are_equal(s, _rth, xi);

        // if we're not near any threshold just divide like normal
        if (!(close_to_sth || close_to_rth)) return ksf_inhomogeneity(i, s)/nus;
        
        double eps = _settings._expansion_eps;
        double s_exp = (close_to_sth) ? _sth : _rth;
        // The only time we expand below (i.e. s_exp - eps) is if s \simeq rth
        bool expand_below = (close_to_rth) && (s <= _rth);

        std::array<complex,3> ecs = rthreshold_expansion(i, s_exp, expand_below);
        return (ecs[0]+ecs[1]*(s-s_exp)+ecs[2]*pow(s-s_exp,2))/nus;
    };

    // Expand the ksf_inhomogeneity near regular thresholds
    std::array<complex,3> raw_iteration::rthreshold_expansion(unsigned int i, double s, bool expand_below)
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

        complex a = (as[0]*f+as[1]*e*fp+as[2]*e*e*fpp)/(as[3]*pow(std::abs(e), n/2.   ));
        complex b = (bs[0]*f+bs[1]*e*fp+bs[2]*e*e*fpp)/(bs[3]*pow(std::abs(e), n/2.+1.));
        complex c = (cs[0]*f+cs[1]*e*fp+cs[2]*e*e*fpp)/(cs[3]*pow(std::abs(e), n/2.+2.));

        if (expand_below) b *= -1;
        return {a,b,c};
    };

    // -----------------------------------------------------------------------
    // The fully regularized integrand also removes the singularity at pth
    // This is only a function of the integration variable (no cauchy kernel)
    complex raw_iteration::regularized_integrand(unsigned int i, double s)
    {
        if (s < _sth)     return NaN<double>();

        // xi is the interval around pth we expand around
        int             n = _n_singularity;
        double         xi = _settings._matching_interval;
        bool close_to_pth = are_equal( s, _pth, xi);
        bool    below_pth = (s <= _pth);

        // Momentum factor we will be dividing out
        complex k = (below_pth) ? sqrt(_pth-s) : +I*sqrt(s-_pth);
        // Expansion coefficients
        std::array<complex,3> ecs = pthreshold_expansion(i, below_pth);

        // if we're not near any threshold just divide like normal
        if (!close_to_pth)
        {
            complex subtracted = half_regularized_integrand(i, s) - half_regularized_integrand(i, _pth);
            // for p-wave also subtract the first derivative
            if (n > 1) subtracted -= (_pth-s)*ecs[0];
            return subtracted/pow(k, n);
        };

        return (below_pth) ? ecs[1]+ecs[2]*sqrt(_pth-s) : I*(ecs[1]+ecs[2]*sqrt(s-_pth));
    };

    // Expand the ksf_inhomogeneity near regular thresholds
    std::array<complex,3> raw_iteration::pthreshold_expansion(unsigned int i, bool expand_below)
    {
        // These coefficients only depend on the order of the singularity
        int n = _n_singularity;
        std::array<int,3> bs, cs, ds;
        switch (n)
        {
            case 1:  // S-waves
            {
                bs = {-3, +3, -2};
                cs = {-3, +3, -2};
                ds = {+3, -4, +4};
                break;
            };
            case 3:  // P-waves
            {
                bs = {-6, +5, -2};
                cs = {+8, -8, +4};
                ds = {+3, -3, +2};
                break;
            };
            default : 
           {
            warning("raw_iteration::pthreshold_expansion", "Invalid singularity order (n="+std::to_string(n)+")!");
            return {NaN<complex>(), NaN<complex>()};
           };
        };

        // Need up to second derivative
        double e = _settings._expansion_eps;
        if (expand_below) e *= -1;
        double s_exp = _pth+e;

        complex fpth = _re_halfreg[i]->Eval(_pth)    + I*_im_halfreg[i]->Eval(_pth);
        complex f    = _re_halfreg[i]->Eval(s_exp)   + I*_im_halfreg[i]->Eval(s_exp);
        complex fp   = _re_halfreg[i]->Deriv(s_exp)  + I*_im_halfreg[i]->Deriv(s_exp);
        complex fpp  = _re_halfreg[i]->Deriv2(s_exp) + I*_im_halfreg[i]->Deriv2(s_exp);

        complex b = (bs[0]*(f-fpth)+bs[1]*e*fp+bs[2]*e*e*fpp)/pow(std::abs(e), (n-1)/2.);
        complex c = (cs[0]*(f-fpth)+cs[1]*e*fp+cs[2]*e*e*fpp)/pow(std::abs(e),  n   /2.);
        complex d = (ds[0]*(f-fpth)+ds[1]*e*fp+ds[2]*e*e*fpp)/pow(std::abs(e), (n+1)/2.);

        if (!expand_below) { c *= -1; };
        return {b, c, d};
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

        // Now we need to evaluate the inhomogenous integral

        // bool no_problem =   (std::real(s) <= _sth - _settings._matching_interval 
        //                   || std::imag(s) >         _settings._infinitesimal);
        // if (no_problem) return polynomial + disperse_with_pth(i, s, {_sth, _cutoff});

        // double p    = (sth + pth)/2;

        // bool below_pth  = (std::real(s) <= _pth);
        // if (below_pth)  return polynomial + disperse_below_pth(i, s);

        // return polynomial + disperse_above_pth(i, s);
        return 0;
    };

    // -----------------------------------------------------------------------

    // Q functions do not have the Cauchy singularity and just need to be regularized in terms of pth
    complex raw_iteration::Q(int n, complex s, std::array<double,2> bounds)
    {
        if (n < 0) return 0.;
        if (!is_odd(n)) return NaN<complex>();

        double x   = bounds[0], y = bounds[1], z = _pth;

        complex argx = (csqrt(z-s)+csqrt(z-x))/(csqrt(z-s)-csqrt(z-x));
        complex argy = csqrt(y-z)/csqrt(z-s);
        complex Q1 = (log(argx)-2*I*atan(argy))/csqrt(z-s);

        if (n == 1) return Q1;

        complex Q3 = (-2*(I/csqrt(y-z)+1/csqrt(z-x)) + Q1)/(z-s);

        if (n == 3) return Q3;

        return error("raw_iteration::Q: Requested Q_n/2 hasnt been added!", NaN<complex>());
    };

        
    // R functions DO have the Cauchy singularity
    complex raw_iteration::R(int n, complex s, std::array<double,2> bounds)
    {
        if (!is_odd(n)) return NaN<complex>();
        
        int pm = (std::imag(s) >= 0) ? +1 : -1;
        double x   = bounds[0], y = bounds[1], z = _kinematics->pth();

        complex argx = (csqrt(z-x)+csqrt(z-s))/(csqrt(z-x)-csqrt(z-s));
        complex argy = (csqrt(z-s)-csqrt(z-y))/(csqrt(z-s)+csqrt(z-y));
        complex R1 = (log(argx)+log(argy)+I*pm*PI)/csqrt(z-s);

        if (n == 1) return R1;

        complex R3 = (2*(1/csqrt(z-y)-1/csqrt(z-x)) + R1)/(z-s);

        if (n == 3) return R3;

        return error("raw_iteration::R: Requested R_n/2 hasnt been added!", NaN<complex>());
    };


    // -----------------------------------------------------------------------
    // Evaluate the integral in different forms depending on where s is

    // If s is away from the singularities, just disperse naively
    complex raw_iteration::disperse_with_pth(unsigned int i, complex s, std::array<double,2> bounds)
    {

        return 0.;

        // int n = _n_singularity;
        // complex N  = _num_pth[i], Np = _dnum_pth[i];
        // auto fdx = [this,i,n,N,Np,s](double x)
        // {
        //     complex subtract = /* (n>1) ? N + Np*(_pth-x) : */ N;
        //     complex integrand;
        //     integrand  = numerator(i,x) - subtract;
        //     integrand /= pow(csqrt(_pth-x+IEPS), n)*(x-s); 
        //     return integrand;
        // };
        // complex integral = boost::math::quadrature::gauss_kronrod<double,61>::integrate(fdx, _sth, _cutoff, 
        //                                                                                 _settings._dispersion_integrator_depth, 1.E-9, NULL);
        // std::array<double,2> bds = {_sth, _cutoff};

        // return integral + N*Q(n,s,bds)/*  + Np*Q(n-2,s,bds) */;
    };

    complex raw_iteration::disperse_with_cauchy(unsigned int i, complex s, std::array<double,2> bounds)
    {
        return 0.;
    };
}; // namespace iterateKT