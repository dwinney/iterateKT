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
    // Initialize by saving relevant quantities and setting up interpoaltions
    void raw_iteration::initialize(basis_grid & dat)
    {
        _sth = _kinematics->sth(); _pth = _kinematics->pth(); _rth = _kinematics->rth();

        // Load up the interpolators from the input data
        for (int i = 0; i < dat.N_basis(); i++)
        {
            _re_inhom.push_back(new ROOT::Math::Interpolator(dat._s_list, dat._re_list[i]));
            _im_inhom.push_back(new ROOT::Math::Interpolator(dat._s_list, dat._im_list[i]));
        };

        // Also half-regularize them and save as an interpolations
        for (int i = 0; i < dat.N_basis(); i++)
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
    // -----------------------------------------------------------------------
    // The discontinuity of the inhomogenous contribution
    // We were fed this from the constructor so we just access the interpolation
    complex raw_iteration::ksf_inhomogeneity(unsigned int i, double s)
    {
        if (s <= _sth || s >= _settings._cutoff || _zeroth) return 0.;
        return _re_inhom[i]->Eval(s) + I*_im_inhom[i]->Eval(s);
    };

    // -----------------------------------------------------------------------
    // Take the above inhomogeneity and divide by the appropriate number of 
    // nu given by n_singularity
    complex raw_iteration::half_regularized_integrand(unsigned int i, double s)
    {
        if (s <= _sth)    return 0.;
        if (_initialized) return _re_halfreg[i]->Eval(s) + I*_im_halfreg[i]->Eval(s);

        // If no interpolation is saved yet, calculate divide by nu explicitly
        int   n = _n_singularity;
        auto xi = _settings._matching_intervals;

        bool close_to_sth = are_equal(s, _sth, xi[0]);
        bool close_to_rth = are_equal(s, _rth, xi[2]);
        complex nus = pow(_kinematics->nu(s, xi), n);

        // if we're not near any threshold just divide like normal
        if (!close_to_sth && !close_to_rth) return ksf_inhomogeneity(i, s)/nus;

        double s_exp = (close_to_sth) ? _sth : _rth;
        double eps   = (close_to_sth) ? _settings._expansion_offsets[0] : _settings._expansion_offsets[2];
        // The only time we expand below (i.e. s_exp - eps) is if s \simeq rth
        if ((close_to_rth) && (s <= _rth)) eps *= -1;

        std::array<complex,3> ecs = rthreshold_expansion(i, s_exp, eps);
        return (ecs[0]+ecs[1]*(s-s_exp)+ecs[2]*pow(s-s_exp,2))/nus;
    };

    // Expand the ksf_inhomogeneity near regular thresholds
    std::array<complex,3> raw_iteration::rthreshold_expansion(unsigned int i, double s, double e)
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

        double s_exp = s+e;
        complex f    = _re_inhom[i]->Eval(s_exp)   + I*_im_inhom[i]->Eval(s_exp);
        complex fp   = _re_inhom[i]->Deriv(s_exp)  + I*_im_inhom[i]->Deriv(s_exp);
        complex fpp  = _re_inhom[i]->Deriv2(s_exp) + I*_im_inhom[i]->Deriv2(s_exp);

        complex a = (as[0]*f+as[1]*e*fp+as[2]*e*e*fpp)/(as[3]*pow(std::abs(e), n/2.   ));
        complex b = (bs[0]*f+bs[1]*e*fp+bs[2]*e*e*fpp)/(bs[3]*pow(std::abs(e), n/2.+1.));
        complex c = (cs[0]*f+cs[1]*e*fp+cs[2]*e*e*fpp)/(cs[3]*pow(std::abs(e), n/2.+2.));

        if (e < 0) b *= -1;
        return {a,b,c};
    };

    // -----------------------------------------------------------------------
    // The fully regularized integrand also removes the singularity at pth
    // This is only a function of the integration variable (no cauchy kernel)
    complex raw_iteration::regularized_integrand(unsigned int i, double s)
    {
        if (s < _sth)    return 0.;

        // xi is the interval around pth we expand around
        int             n = _n_singularity;
        double         xi = _settings._matching_intervals[1];
        bool close_to_pth = are_equal( s, _pth, xi);
        bool    below_pth = (s < _pth);

        // Momentum factor we will be dividing out
        double eps = (!below_pth - below_pth)*_settings._expansion_offsets[1];
        complex k  = (below_pth) ? csqrt(_pth-s) : I*csqrt(s-_pth);

        // Expansion coefficients
        std::array<complex,4> ecs = pthreshold_expansion(i, eps);

        // if we're not near any threshold just divide like normal
        if (!close_to_pth)
        {
            complex subtracted = half_regularized_integrand(i, s) - ecs[0];
            // for p-wave also subtract the first derivative
            if (n > 1) subtracted -= (_pth-s)*ecs[1];

            return subtracted/pow(k, n);
        };

        return ecs[2]+ecs[3]*k;
    };

    // Expand the ksf_inhomogeneity near regular thresholds
    std::array<complex,4> raw_iteration::pthreshold_expansion(unsigned int i, double epsilon)
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

        if (epsilon < 0)
        {
            bs[0] *= -1; bs[2] *= -1;
            cs[0] *= -1; cs[2] *= -1;
            ds[1] *= -1;
        };

        // Need up to second derivative
        double s_exp = _pth+epsilon;

        double e = abs(epsilon);
        complex a    = _re_halfreg[i]->Eval(_pth)    + I*_im_halfreg[i]->Eval(_pth);
        complex f    = _re_halfreg[i]->Eval(s_exp)   + I*_im_halfreg[i]->Eval(s_exp);
        complex fp   = _re_halfreg[i]->Deriv(s_exp)  + I*_im_halfreg[i]->Deriv(s_exp);
        complex fpp  = _re_halfreg[i]->Deriv2(s_exp) + I*_im_halfreg[i]->Deriv2(s_exp);

        complex b = (bs[0]*(f-a)+bs[1]*e*fp+bs[2]*e*e*fpp)/pow(e, (n-1.)/2.);
        complex c = (cs[0]*(f-a)+cs[1]*e*fp+cs[2]*e*e*fpp)/pow(e,  n   /2.);
        complex d = (ds[0]*(f-a)+ds[1]*e*fp+ds[2]*e*e*fpp)/pow(e, (n+1.)/2.);
        
        if (epsilon > 0) c *= -I;
        return {a, b, c, d};
    };

    // -----------------------------------------------------------------------
    // Evaluate the `basis' function. This deviates from the typical 
    // defintion by not including the overall factor of the omnes function

    complex raw_iteration::integral(unsigned int i, complex sc)
    {
        if (is_zero(sc) || _zeroth) return 0.;
        
        // If we're sufficiently far from pth we can just integrate without issue
        bool no_problem = (std::real(sc) < _sth || abs(std::imag(sc)) > _settings._infinitesimal);
        if (no_problem) return disperse_with_pth(i, sc, {_sth, _settings._cutoff});

        // If we're too close to the real line, we evalaute with ieps perscriptions
        complex ieps = sign(std::imag(sc)) * I*_settings._infinitesimal;
        double s     = std::real(sc);

        // Split the integrals into two pieces, one which contains the pth singularity
        // and another with the cauchy singularity
        double p    = (s + _pth)/2;
        std::array<double,2> lower = {_sth, p}, upper = {p, _settings._cutoff};

        // If we are evaluating exactly at this point we have a problem
        // If we're too close to pth just evaluate above and below it and linear interpolate
        double xi = _settings._matching_intervals[1];
        bool on_pth = are_equal(s, _pth, xi/2. * 0.9);
        if (on_pth)
        {
            double  x1 = _pth - xi/2., x2 = _pth + xi/2.;
            complex f1 = integral(i, x1+ieps), f2 = integral(i, x2+ieps);
            return  f1 + (f2 - f1)*(s - x1)/(x2 - x1);
        };

        // Else evaluate the integrals
        bool below_pth  = (p <= _pth);
        if (below_pth) return disperse_with_cauchy(i, s+ieps, lower)+ disperse_with_pth(i, s+ieps, upper);        
        return disperse_with_pth(i, s+ieps, lower) + disperse_with_cauchy(i, s+ieps, upper);
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
    complex raw_iteration::R(int n, complex sc, std::array<double,2> bounds)
    {
        if (!is_odd(n)) return NaN<complex>();
        
        int pm   = (std::imag(sc) <= 0) ? +1 : -1;
        double s = std::real(sc);
        if (are_equal(s,_sth)) return 0.;
        
        double x   = bounds[0], y = bounds[1], z = _kinematics->pth();

        complex argx = (csqrt(z-x)+csqrt(z-s))/(csqrt(z-x)-csqrt(z-s));
        complex argy = (csqrt(z-s)-csqrt(z-y))/(csqrt(z-s)+csqrt(z-y));
        complex R1 = (log(argx)+log(argy)+I*pm*PI)/csqrt(z-s);

        if (n == 1) return R1;

        complex R3 = (2*(1/csqrt(z-y)-1/csqrt(z-x)) + R1)/std::abs(z-s);

        if (n == 3) return R3;

        return error("raw_iteration::R: Requested R_n/2 hasnt been added!", NaN<complex>());
    };


    // -----------------------------------------------------------------------
    // Evaluate the integral in different forms depending on where s is

    // This is the integral which contains and handles the pth singularity
    complex raw_iteration::disperse_with_pth(unsigned int i, complex s, std::array<double,2> bounds)
    {
        using namespace boost::math::quadrature;

        // We'll need to calculate the 
        int       n = _n_singularity;
        int      pm = (std::real(s) < _pth) ? -1 : +1;
        auto coeffs = pthreshold_expansion(i, pm*_settings._expansion_offsets[1]);
        complex a = coeffs[0], b = coeffs[1];
        
        auto fdx = [this,i,s](double x){ return regularized_integrand(i,x)/(x-s); };
        complex integral = (_settings._adaptive_dispersion) ? gauss_kronrod<double,61>::integrate(fdx, bounds[0], bounds[1], _settings._dispersion_depth, 1.E-9, NULL)
                                                            : gauss<double,N_GAUSS>::   integrate(fdx, bounds[0], bounds[1]);

        return integral + a*Q(n,s,bounds) + b*Q(n-2,s,bounds);
    };

    // This is the integral which contains and handles the cauchy singularity
    complex raw_iteration::disperse_with_cauchy(unsigned int i, complex s, std::array<double,2> bounds)
    {
        using namespace boost::math::quadrature;

        int n = _n_singularity;
        complex a = half_regularized_integrand(i, std::real(s));
        auto fdx = [this,i,s,a,n](double x)
        { 
            complex kx = (x <= _pth) ? csqrt(_pth - x) : I*csqrt(x - _pth);
            return (half_regularized_integrand(i,x) - a)/(x-s)/pow(kx,n); 
        };
        complex integral = (_settings._adaptive_dispersion) ? gauss_kronrod<double,61>::integrate(fdx, bounds[0], bounds[1], _settings._dispersion_depth, 1.E-9, NULL)
                                                            : gauss<double,N_GAUSS>::   integrate(fdx, bounds[0], bounds[1]);

        return integral + a*R(n,s,bounds);
    };
}; // namespace iterateKT