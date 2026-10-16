
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
    // Set up the zeroth iteration 
    // and set up the s values which will be used for interpolations later on
    void raw_isobar::initialize()
    {
        // Add the 'zeroth' iteration to the list
        _iterations.push_back(new_iteration());
        _ieps =  I*_settings._infinitesimal;
        double sth = _kinematics->sth(), pth = _kinematics->pth();

        // We interpolate different regions with different graining
        // Each s[i] marks a boundary of a region interpolated with N[i] points
        // eps is added to make a gap between regions and ensure each
        // interval is monotonic increasing
        double eps   = _settings._interpolation_offset;
        double xi    = _settings._expansion_offsets[1];
        double s_cut = _settings._cutoff;
        double s_int = _settings._intermediate_energy;

        int order = (pth > s_int) + (pth > s_cut);
        std::array<double, 5> s = {sth, pth-1.5*xi, pth+1.5*xi, s_cut, s_int};
        std::sort(s.begin(), s.end());

        std::array<int,3> N = _settings._interpolation_points;
        std::array<int,4> n;
        switch (order)
        {
            case 0: // {sth, pth-, pth+, s_int, s_cut}
            {
                n[0] = int((s[1]-s[0])/(s[3]-s[0])*N[0]); // s0 -> s1
                n[1] = N[1];                              // s1 -> s2
                n[2] = N[0] - n[0];                       // s2 -> s3
                n[3] = N[2];                              // s3 -> s4
                break;
            };
            case 1: // {sth, s_int, pth-, pth+, s_cut}
            {
                n[0] = N[0];
                n[1] = int((s[2]-s[1])/(s[4]-s[1])*N[2]);
                n[2] = N[1];
                n[3] = N[2] - n[1];
                break;
            };
            case 2: // {sth, s_int, s_cut, pth-, pth+}
            {
                n[0] = N[0];
                n[1] = N[2];
                n[2] = 0;
                n[3] = 0;
                break;
            };
            default: fatal("isobar::initialize()", "Weird error! Don't know how you got here");
        }

        // Populate the s values we need to evaluate at
        for (int i = 0; i < 4; i++)
        {
            double s1 = s[i]+(i!=0)*eps, s2 = s[i+1];
            if (s1 >= s2 || n[i] < 2) continue;
            for (int j = 0; j < n[i]; j++)
            {
                double x = s1 + j*(s2-s1)/(n[i]-1);
                _s_list.push_back(x);
            }
        };

        // If we cutoff before pth we dont need to exclude it
        if (order == 2) return;

        // We also need to be able to exclude a part of the isobars around pth
        int n_ex = _settings._exclusion_points/2;
        double low, high;
        low  = pth - 2*_settings._exclusion_offsets[0];
        high = pth -   _settings._exclusion_offsets[0];
        for (int k = 0; k <= n_ex; k++) _s_around_pth.push_back(low-k*(low-high)/n_ex);
        low  = pth +   _settings._exclusion_offsets[1];
        high = pth + 2*_settings._exclusion_offsets[1];
        for (int k = 0; k <= n_ex; k++) _s_around_pth.push_back(low-k*(low-high)/n_ex);
    };

    // Save an interpolation of the LHC since this never changes and is called a lot
    void raw_isobar::interpolate_lhc()
    {
        std::vector<double> lhc;
        for (auto s : _s_list) lhc.push_back( sin(this->phase_shift(s))/std::abs(omnes(s+_ieps)) );
        _lhc.SetData(_s_list, lhc); 
        _lhc_interpolated = true;
    };

    // ----------------------------------------------------------------------- 
    // Evaluate the Omnes functions from the given phase_shift

    complex raw_isobar::omnes(complex s)
    {
        using namespace boost::math::quadrature;

        // Only explciitly evaluate above cut, below is given by Schwarz
        if (imag(s) < 0) return conj(omnes(conj(s)));

        // bounds of integration
        double low  = _kinematics->sth();
        double mid  = _settings._intermediate_energy;
        double high = _settings._omnes_cutoff;
        
        // If we're sufficiently far from the real axis just do the integral naively
        if (imag(s) > imag(_ieps))
        {
            auto fdx = [this, s](double x)
            {
                complex integrand;
                integrand  = phase_shift(x);
                integrand *= (s/x); // One subtraction
                integrand /= (x-s); 
                return integrand;
            };
            // If using gauss gauss-legendre, split the integral into two pieces to avoid systematic errors at low energies
            complex integral = gauss_kronrod<double,N_GAUSS_OMNES>::integrate(fdx, low, mid,  _settings._omnes_integrator_depth, 1.E-9, NULL) 
                             + gauss_kronrod<double,N_GAUSS_OMNES>::integrate(fdx, mid, high, _settings._omnes_integrator_depth, 1.E-9, NULL);
            return exp(integral/M_PI);
        }

        // If we're close to the real axis, we split the integration in two parts
        // to properly handle the Principle Value and ieps perscription

        double RHCs = (real(s) <= _kinematics->sth()) ? 0 : phase_shift(real(s));
        auto    fdx = [this,s,RHCs](double x)
        {
            complex integrand;
            integrand  = phase_shift(x) - RHCs;
            integrand *= (s/x); // One subtraction
            integrand /= (x-(s+_ieps)); 
            return integrand;
        };

        // If using gauss gauss-legendre, split the integral into two pieces to avoid systematic errors at low energies
        complex integral = gauss_kronrod<double,N_GAUSS_OMNES>::integrate(fdx, low, mid,  _settings._omnes_integrator_depth, 1.E-9, NULL) 
                         + gauss_kronrod<double,N_GAUSS_OMNES>::integrate(fdx, mid, high, _settings._omnes_integrator_depth, 1.E-9, NULL);

        complex logarithm = RHCs * log(1.-(s+_ieps) / low);
        return exp((integral-logarithm)/M_PI);
    };

    // ----------------------------------------------------------------------- 
    // Basis funciton is given by the form
    // Omega(s)*(P(s) + I(s))

    // Specify a given iteration to use when outputting the basis_function
    complex raw_isobar::basis_function(unsigned int iter_id, unsigned int basis_id, complex s)
    {
        if (iter_id  >=  _iterations.size())       return 0;
        if (basis_id >= _subtractions->N_basis())  return 0;

        bool no_P = (_subtractions->get_id(basis_id) != get_id());
        complex P = (no_P) ? 0 : _subtractions->driving_term(basis_id, s);
        if ( is_zero(s) ) return P;
        return omnes(s)*(P + pow(s,_max_sub)/PI*_iterations[iter_id]->integral(basis_id, s));
    };

    // Without an iter_id we just take the latest iteration
    complex raw_isobar::basis_function(unsigned int basis_id, complex x)
    { 
        return basis_function(_iterations.size()-1, basis_id, x); 
    };
    
    // ----------------------------------------------------------------------- 
    // Take the saved interpolation settings and output the necessary arrays
    
    basis_grid raw_isobar::calculate_next(std::vector<isobar> & previous)
    {
        basis_grid output;
        output._n_singularity = 2*angular_momentum()+1;
        output._s_list        = _s_list;
        output._s_around_pth  = _s_around_pth;
        
        double x = _settings._iteration_rate_intercept + _iterations.size()*_settings._iteration_rate_slope;
        if (x >= 1) x = 1;

        // Sum over basis functions
        for (int i = 0; i < _subtractions->N_basis(); i++)
        {
            std::vector<double> re, im;
            for (auto s : _s_list)
            {
                complex new_disc = LHC(s)/pow(s,_max_sub)*pinocchio_integral(i,s,previous);
                complex old_disc = _iterations.back()->ksf_inhomogeneity(i, s);
                complex weighted = x*new_disc - (1-x)*old_disc;
                re.push_back( real(weighted) );
                im.push_back( imag(weighted) );
            };
            output._re_list.push_back(re);
            output._im_list.push_back(im);
        };

        return output;
    };

    // Filter which region of the pinnochio we are evalutating at and call the appropriate function
    complex raw_isobar::pinocchio_integral(unsigned int basis_id, double s, std::vector<isobar> & previous)
    {
        int region = (s > _kinematics->B()) + (s > _kinematics->C()) + (s > _kinematics->D());
        switch (region)
        {
            // Both s+ and s- real and above cut
            case 0: 
            case 3:
            {
                double sp = real(_kinematics->t_plus(s));
                double sm = real(_kinematics->t_minus(s));

                // If region 0, we needs an ieps to avoid cuts
                // in region 3, both sp & sm are negative and dont need
                double ieps = (region == 0) ? +1 : 0;

                // Check if we have cusp
                double sc = _settings._extra_cusp;
                if (sp > sc && sm < sc)
                {
                    return linear_segment(basis_id, {sm, sc, ieps}, s, previous) 
                         + linear_segment(basis_id, {sc, sp, ieps}, s, previous);
                };
                return linear_segment(basis_id, {sm, sp, ieps}, s, previous);
            };
            // s+ is above cut but s- is below cut
            case 1:
            {
                complex integ_above, integ_below;
                double sth = _kinematics->sth(), sc = _settings._extra_cusp;

                // Check if we have cusp
                double sp = real(_kinematics->t_plus(s));
                if (sp > sc && sc > sth)
                {
                    integ_above = linear_segment(basis_id, {sth, sc, +1}, s, previous) 
                                + linear_segment(basis_id, {sc,  sp, +1}, s, previous);
                }
                else integ_above = linear_segment(basis_id, {sth, sp, +1}, s, previous);

                // Do same with the segment below cut
                double sm = real(_kinematics->t_minus(s));
                if (sm > sc && sc > sth)
                {
                    integ_below = linear_segment(basis_id, {sm, sc,   -1}, s, previous) 
                                + linear_segment(basis_id, {sc,  sth, -1}, s, previous); 
                }
                else integ_below = linear_segment(basis_id, {sm, sth, -1}, s, previous);

                return integ_above + integ_below;
            };
            // In the curved "egg" portion
            case 2:  return curved_segment(basis_id, s, previous);
            default: return 0.;
        };

        return NaN<complex>();
    };

    // Integrate along a linear segment +ieps above the real axis
    complex raw_isobar::linear_segment(unsigned int basis_id, std::array<double,3> bounds, double s, std::vector<isobar> & previous_list)
    {
        using namespace boost::math::quadrature;

        //  sum over all the previous isobars with their appropriate kernels
        double pm = bounds[2];
        auto fdx = [this,previous_list,s,basis_id,pm](double t)
        {
            complex sum = 0;
            for (auto previous : previous_list) 
            {
                complex K = ksf_kernel(previous->get_id(),s,t);
                if (is_zero(K)) continue;
                sum += K*previous->basis_function(basis_id,t+pm*_ieps);
            };
            return sum;
        };
        return gauss_kronrod<double,N_GAUSS_ANGULAR>::integrate(fdx, bounds[0], bounds[1], _settings._angular_integrator_depth, 1.E-9, NULL);
    };

    // Integrate along the curved segment of pinnochio's head
    complex raw_isobar::curved_segment(unsigned int basis_id, double s, std::vector<isobar> & previous_list)
    {
        using namespace boost::math::quadrature;

        //  sum over all the previous isobars with their appropriate kernels
        complex ieps = I*_settings._infinitesimal;
        auto fdx = [this,previous_list,s,ieps,basis_id](double phi)
        {
            complex sum = 0;
            for (auto previous : previous_list) 
            {
                complex t = this->_kinematics->t_curve(phi);
                complex K = ksf_kernel(previous->get_id(),s,t);
                if (is_zero(K)) continue;
                sum += K*previous->basis_function(basis_id,t);
            };
            sum *= this->_kinematics->jacobian(phi);
            return sum;
        };
        return gauss_kronrod<double,N_GAUSS_ANGULAR>::integrate(fdx, _kinematics->phi_minus(s), _kinematics->phi_plus(s), _settings._angular_integrator_depth, 1.E-9, NULL);
    };
}; // namespace iterateKT