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
        using ROOT::Math::Interpolator;

        _sth = _kinematics->sth(); _pth = _kinematics->pth(); _rth = _kinematics->rth();

        auto interp_type = _settings._interpolation_type;
        // Load up the interpolators from the input data
        for (int i = 0; i < dat.N_basis(); i++)
        {
            _re.push_back(new Interpolator(dat._s_list, dat._re_list[i], interp_type));
            _im.push_back(new Interpolator(dat._s_list, dat._im_list[i], interp_type));

            // Also calculate all the expansion coefficients so we dont need to calculate
            // derivatives multiple times
            double sth_eps = _settings._expansion_offsets[0];
            _sth_expansion.push_back(rthreshold_expansion(i, _sth, +sth_eps));

            double pth_eps = _settings._expansion_offsets[1];
            _below_pth_expansion.push_back(pthreshold_expansion(i, -pth_eps));
            _above_pth_expansion.push_back(pthreshold_expansion(i, +pth_eps));

            double rth_eps = _settings._expansion_offsets[2];
            _below_rth_expansion.push_back(rthreshold_expansion(i, _rth, -rth_eps));
            _above_rth_expansion.push_back(rthreshold_expansion(i, _rth, +rth_eps));
        };

        // Finally, because we have singularities at pth, evaluate below and above it
        // and interpolate in between to get a smooth curve
        for (int i = 0; i < dat.N_basis(); i++)
        {
            std::vector<double> rep, imp, rem, imm;
            for (auto s : dat._s_around_pth)
            {
                complex disp_plus = integral(i, s+IEPS);
                rep.push_back( real(disp_plus) );
                imp.push_back( imag(disp_plus) );

                complex disp_minus = integral(i, s-IEPS);
                rem.push_back( real(disp_minus) );
                imm.push_back( imag(disp_minus) );
            };
            // Use AKIMA interpolation instead of cspline
            auto exc_interp_type = ROOT::Math::Interpolation::Type::kAKIMA;
            _re_excluded_above.push_back(new Interpolator(dat._s_around_pth, rep, exc_interp_type));
            _im_excluded_above.push_back(new Interpolator(dat._s_around_pth, imp, exc_interp_type));
            _re_excluded_below.push_back(new Interpolator(dat._s_around_pth, rem, exc_interp_type));
            _im_excluded_below.push_back(new Interpolator(dat._s_around_pth, imm, exc_interp_type));
        };

        _initialized = true;
    };
    // -----------------------------------------------------------------------
    // The discontinuity of the inhomogenous contribution
    // We were fed this from the constructor so we just access the interpolation
    complex raw_iteration::ksf_inhomogeneity(unsigned int i, double s)
    {
        if (i >= _re_inhom.size()) fatal("raw_iteration", "Requested out of scope basis function!");
        if (s <= _sth || s >= _settings._cutoff || _zeroth) return 0.;
        return _re[i]->Eval(s) + I*_im[i]->Eval(s);
    };

    // -----------------------------------------------------------------------
    // Take the above inhomogeneity and divide by the appropriate number of 
    // nu given by n_singularity
    complex raw_iteration::half_regularized_integrand(unsigned int i, double s)
    {
        if (are_equal(s, _sth, _settings._infinitesimal)) return 0.;

        auto     xi = _settings._matching_intervals;
        complex nus = pow(_kinematics->nu(s, xi), _n);
        
        bool close_to_sth = are_equal(s, _sth, xi[0]);
        bool close_to_rth = are_equal(s, _rth, xi[2]);
        // if we're not near any threshold just divide like normal
        if (!close_to_sth && !close_to_rth) return ksf_inhomogeneity(i, s)/nus;

        double s_exp; std::array<complex,3> coeffs;
        if (close_to_sth){ s_exp = _sth; coeffs = _sth_expansion[i]; }
        else             { s_exp = _rth; coeffs = (s > _rth) ? _above_rth_expansion[i] : _below_rth_expansion[i]; };

        complex expansion = 0;
        for (int i = 0; i < coeffs.size(); i++) expansion += coeffs[i] * pow(abs(s-s_exp), i);
        return expansion / nus;
    };

    // Expand the ksf_inhomogeneity near regular thresholds
    std::array<complex,3> raw_iteration::rthreshold_expansion(unsigned int i, double s, double e)
    {
        // These coefficients only depend on the order of the singularity
        std::array<int,4> as, bs, cs;
        switch (_n)
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
            case 7: // F-waves
            {
                as = {+99, -36, +4, 8};
                bs = {-77, +32, -4, 4};
                cs = {+63, -28, +4, 8};
                break;
            };
            default : 
           {
            warning("raw_iteration::expansion_coefficients", "Invalid singularity order (n="+to_string(_n)+")!");
            return {NaN<complex>(), NaN<complex>(), NaN<complex>()};
           };
        };

        double s_exp = s+e;
        complex f    = _re[i]->Eval(s_exp)   + I*_im[i]->Eval(s_exp);
        complex fp   = _re[i]->Deriv(s_exp)  + I*_im[i]->Deriv(s_exp);
        complex fpp  = _re[i]->Deriv2(s_exp) + I*_im[i]->Deriv2(s_exp);

        complex a = (as[0]*f+as[1]*e*fp+as[2]*e*e*fpp)/(as[3]*pow(abs(e), _n/2.   ));
        complex b = (bs[0]*f+bs[1]*e*fp+bs[2]*e*e*fpp)/(bs[3]*pow(abs(e), _n/2.+1.));
        complex c = (cs[0]*f+cs[1]*e*fp+cs[2]*e*e*fpp)/(cs[3]*pow(abs(e), _n/2.+2.));

        return {a, b, c};
    };

    // -----------------------------------------------------------------------
    // The fully regularized integrand also removes the singularity at pth
    // This is only a function of the integration variable (no cauchy kernel)
    complex raw_iteration::regularized_integrand(unsigned int i, double s)
    {
        if (s < _sth)    return 0.;

        // xi is the interval around pth we expand around
        double         xi = _settings._matching_intervals[1];
        bool close_to_pth = are_equal( s, _pth, xi);

        // Sign of the s_exp = s \pm epsilon expansion
        auto coeffs =  (s < _pth) ? _below_pth_expansion[i] : _above_pth_expansion[i];

        // Momentum factor we will be dividing out
        complex ks = k(s);

        // if we're not near any threshold just divide like normal
        if (!close_to_pth)
        {
            complex subtracted = half_regularized_integrand(i, s);
            // Subtract away as many terms as required to cancel the singularity at pth
            for (int i = 0; i <= _l; i++) subtracted -= coeffs[i] * pow(_pth-s,i);
            return subtracted/pow(ks, _n);
        };

        // Finally, if we're close to pth return the expansion in k(s)
        complex expansion = 0;
        for (int i = 1+_l; i < coeffs.size(); i++) expansion += coeffs[i] * pow(ks, i-1-_l);
        return expansion;
    };

    // Expand the ksf_inhomogeneity near regular thresholds
    std::array<complex,4> raw_iteration::pthreshold_expansion(unsigned int i, double epsilon)
    {
        // These coefficients only depend on the order of the singularity
        // and whether we're above or below pth
        std::array<int,3> bs, cs, ds;
        switch (sign(epsilon)*int(_n))
        {
            case +1:  // S-waves (above pth)
            {
                bs = {-3, +3, -2};
                cs = {+3, -4, +4};
                ds = {+1, -1, +2};
                break;
            };
            case -1:  // S-waves (below pth)
            {
                bs = {+3, +3, +2};
                cs = {-3, -4, -4};
                ds = {+1, +1, +2};
                break;
            };
            case +3:  // P-waves (above pth)
            {
                bs = {-6, +5, -2};
                cs = {-8, +8, -4};
                ds = {+3, -3, +2};
                break;
            };
            case -3:  // P-waves (below pth)
            {
                bs = {+6, +5, +2};
                cs = {-8, -8, -4};
                ds = {+3, +3, +2};
                break;
            };
            default : fatal("raw_iteration::pthreshold_expansion", "Invalid singularity order (n="+std::to_string(_n)+")!");
        };
        
        // Need up to second derivative
        auto    F     = [this,i](double s){ return half_regularized_integrand(i,s); };
        complex f0    = F(_pth);
        complex f     = F(_pth+epsilon);
        complex fp    = central_difference_derivative<complex>(1, F, _pth+epsilon, _settings._derivative_h);
        complex fpp   = central_difference_derivative<complex>(2, F, _pth+epsilon, _settings._derivative_h);
       
        double  e = abs(epsilon);
        complex a = f0;
        complex b = (bs[0]*(f-f0)+bs[1]*e*fp+bs[2]*e*e*fpp)/pow(e, (_l+1)/2.);
        complex c = (cs[0]*(f-f0)+cs[1]*e*fp+cs[2]*e*e*fpp)/pow(e, (_l+2)/2.);
        complex d = (ds[0]*(f-f0)+ds[1]*e*fp+ds[2]*e*e*fpp)/pow(e, (_l+3)/2.);

        return {a, b, c, d};
    };

    // -----------------------------------------------------------------------
    // Evaluate the `basis' function. This deviates from the typical 
    // defintion by not including the overall factor of the omnes function

    complex raw_iteration::integral(unsigned int i, complex sc)
    {
        if (is_zero(sc) || _zeroth) return 0.;
        
        // If we're sufficiently far from pth we can just integrate without issue
        bool no_problem = (real(sc) < _sth || abs(imag(sc)) > _settings._infinitesimal);
        if  (no_problem) return disperse_with_pth(i, sc, {_sth, _settings._cutoff});

        // If we're too close to the real line, we evalaute with ieps perscriptions
        double s = real(sc), eps = sign(imag(sc))*_settings._infinitesimal;

        // Split the integrals into two pieces, one which contains the pth singularity
        // and another with the cauchy singularity
        double p    = (s + _pth)/2.;
        std::array<double,2> lower = {_sth, p}, upper  = {p, _settings._cutoff};

        // If we are evaluating exactly at this point we have a problem
        // If we're too close to pth just evaluate above and below it and linear interpolate
        bool in_excluded = (s >= _pth - _settings._exclusion_offsets[0]) &&
                           (s <= _pth + _settings._exclusion_offsets[1]) && _initialized;
        if (in_excluded) return (eps > 0) ? _re_excluded_above[i]->Eval(s) + I*_im_excluded_above[i]->Eval(s)
                                          : _re_excluded_below[i]->Eval(s) + I*_im_excluded_below[i]->Eval(s);

        // Else evaluate the integrals
        if (p <= _pth) return disperse_with_cauchy(i, s+I*eps, lower) + disperse_with_pth   (i, s+I*eps, upper);
                       return disperse_with_pth   (i, s+I*eps, lower) + disperse_with_cauchy(i, s+I*eps, upper);
    };

    // -----------------------------------------------------------------------

    // Q functions analytically do integrals. 
    // These are of the form \int dx 1\frac{1}{(pth-x)^{n/2}(x-s-ieps)}
    // and can be defined recoursively but recursion is slower than nonrecursion
    complex raw_iteration::Q(int n, complex s, std::array<double,2> bounds)
    {
        if (n <= 0 || n%2==0) return 0.;
        complex ks = k(s), kx = k(bounds[0]), ky = k(bounds[1]);
        if (n == 1) return 2*(atanh(kx/ks)-atanh(ky/ks))/ks;
        return (2/ky-2/kx+Q(n-2,s,bounds))/(_pth-s)/(n-2);
    };

    // -----------------------------------------------------------------------
    // Evaluate the integral in different forms depending on where s is

    // This is the integral which contains and handles the pth singularity
    complex raw_iteration::disperse_with_pth(unsigned int i, complex s, std::array<double,2> bounds)
    {
        using namespace boost::math::quadrature;
        auto fdx = [this,i,s](double x){ return regularized_integrand(i,x)/(x-s); };
        // Integrate on either side of the pth singularity 
        complex integral = gauss_kronrod<double,N_GAUSS_PSEUDO>::integrate(fdx, bounds[0], _pth, _settings._pseudo_integrator_depth, 1.E-9, NULL)
                         + gauss_kronrod<double,N_GAUSS_PSEUDO>::integrate(fdx, _pth, bounds[1], _settings._pseudo_integrator_depth, 1.E-9, NULL);
        // Add back the analytic pieces we subtracted before
        auto coeffs = (real(s) < _pth) ? _below_pth_expansion[i] : _above_pth_expansion[i];
        for (int i = 0; i <= _l; i++) integral += coeffs[i] * Q(_n-2*i,s,bounds);
        return integral;
    };

    // This is the integral which contains and handles the cauchy singularity
    complex raw_iteration::disperse_with_cauchy(unsigned int i, complex s, std::array<double,2> bounds)
    {
        using namespace boost::math::quadrature;
        complex a = half_regularized_integrand(i, real(s));
        auto fdx = [this,i,s,a](double x)
        { 
            return (half_regularized_integrand(i,x) - a)/(x-s)/pow(k(x),_n); 
        };
        complex integral = gauss_kronrod<double,N_GAUSS_CAUCHY>::integrate(fdx, bounds[0], real(s), _settings._cauchy_integrator_depth, 1.E-9, NULL)
                         + gauss_kronrod<double,N_GAUSS_CAUCHY>::integrate(fdx, real(s), bounds[1], _settings._cauchy_integrator_depth, 1.E-9, NULL);
        return integral + a*Q(_n,s,bounds);
    };
}; // namespace iterateKT