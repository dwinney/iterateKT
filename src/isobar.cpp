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
    // Set up the zeroth iteration but also save interpoaltions of the LHC
    // amplitude since this is always the same
    void raw_isobar::initialize()
    {
        // Add the 'zeroth' iteration to the list
        _iterations.push_back(new_iteration(_max_sub));
        _ieps =  I*_settings._infinitesimal;

        // We interpolate different regions with different graining
        // Each si marks a boundary of a region 
        // eps is added to make a gap between regions and ensure each
        // interval is monotonic increasing

        double eps = _settings._interpolation_offset;
        double s0  = _kinematics->sth();
        double s1  = _kinematics->pth() - 1.5*_settings._expansion_offsets[1] - eps; // This region is around the point we will 
        double s2  = _kinematics->pth() + 1.5*_settings._expansion_offsets[1];       // expand the angular average around pth
        double s3  = _settings._intermediate_energy;
        double s4  = _settings._cutoff;

        std::array<int,3> Ns = _settings._interpolation_points;
        int N_0  = Ns[0];
        int N_1  = int((s1 - s0)/(s3 - s0)*N_0);
        int N_2  = Ns[1];
        int N_3  = N_0 - N_1;
        int N_4  = Ns[2];

        for (int i = 0; i < N_1; i++) _s_list.push_back(s0 + (s1 - s0)*double(i)/double(N_1-1));
        s1 += eps;
        for (int i = 0; i < N_2; i++) _s_list.push_back(s1 + (s2 - s1)*double(i)/double(N_2-1));
        s2 += eps;
        for (int i = 0; i < N_3; i++) _s_list.push_back(s2 + (s3 - s2)*double(i)/double(N_3-1));
        s3 += eps;
        for (int i = 0; i < N_4; i++) _s_list.push_back(s3 + (s4 - s3)*double(i)/double(N_4-1));
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
        if (imag(s) < 0) return std::conj(omnes(std::conj(s)));

        // bounds of integration
        double low = _kinematics->sth();
        double high = std::numeric_limits<double>::infinity();
        
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
            complex integral = (_settings._adaptive_omnes) ? gauss_kronrod<double,61>::integrate(fdx, low, high, _settings._omnes_depth, 1.E-9, NULL)
                                                           : gauss<double,N_GAUSS>::   integrate(fdx, low, high);
            return exp(integral/M_PI);
        }

        // If we're close to the real axis, we split the integration in two parts
        // to properly handle the Principle Value and ieps perscription

        double  RHCs = phase_shift(real(s));
        auto fdx = [this,s,RHCs](double x)
        {
            complex integrand;
            integrand  = phase_shift(x) - RHCs;
            integrand *= (s/x); // One subtraction
            integrand /= (x-(s+_ieps)); 
            return integrand;
        };

        complex integral = (_settings._adaptive_omnes) ? gauss_kronrod<double,61>::integrate(fdx, low, high, _settings._omnes_depth, 1.E-9, NULL)
                                                       : gauss<double,N_GAUSS>::   integrate(fdx, low, high);

        complex logarithm = (std::real(s) <= _kinematics->sth()) ? 0 : RHCs * log(1.-(s+_ieps) / low);
    
        return exp((integral-logarithm)/M_PI);
    };

    // ----------------------------------------------------------------------- 
    // Specify a given iteration to use when outputting the basis_function
    complex raw_isobar::basis_function(unsigned int iter_id, unsigned int basis_id, complex x)
    {
        if (iter_id  > _iterations.size()) return error("Requested iteration does not exist!", NaN<complex>());
        if (basis_id > _max_sub - 1)  return error("Requested basis function does not exist!", NaN<complex>());

        return omnes(x)*_iterations[iter_id]->basis_factor(basis_id, x);
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
    
    basis_grid raw_isobar::calculate_next(std::vector<isobar> & previous)
    {
        basis_grid output;
        output._s_list = _s_list;
        // Sum over basis functions
        for (int i = 0; i < _max_sub; i++)
        {
            std::vector<double> re, im;
            for (auto s : _s_list)
            {
                complex ksf_disc = LHC(s)*pinocchio_integral(i,s,previous)/pow(s,_max_sub);
                re.push_back( std::real(ksf_disc) );
                im.push_back( std::imag(ksf_disc) );
            };
            output._re_list.push_back(re);
            output._im_list.push_back(im);
        };

        return output;
    };

    // Filter which region of the pinnochio we are evalutating at and call the appropriate function
    complex raw_isobar::pinocchio_integral(unsigned int basis_id, double s, std::vector<isobar> & previous)
    {
        if (s < _kinematics->A()) return error("isobar::angular_integral", 
                                                "Trying to evaluate angular integral below threshold!", NaN<complex>());
        
        int region = (s > _kinematics->B()) + (s > _kinematics->C()) + (s > _kinematics->D());
        
        switch (region)
        {
            // Both s+ and s- real and above cut
            case 0: 
            case 3:
            {
                double sp = std::real(_kinematics->t_plus(s));
                double sm = std::real(_kinematics->t_minus(s));
                return linear_segment(basis_id, {sm, sp, 0.}, s, previous);
            };
            // s+ is above cut but s- is below cut
            case 1:
            {
                double sp = std::real(_kinematics->t_plus(s));
                double sm = std::real(_kinematics->t_minus(s));
                return linear_segment(basis_id, {_kinematics->sth(), sp, +1}, s, previous)  //+ieps
                     + linear_segment(basis_id, {sm, _kinematics->sth(), -1}, s, previous); //-ieps
            };
            // In the curved "egg" portion
            case 2: return curved_segment(basis_id, s, previous);
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
            sum+=kernel(previous->id(),s,t)*previous->basis_function(basis_id,t+pm*_ieps);
            return sum;
        };
        return (_settings._adaptive_angular) ? gauss_kronrod<double,61>::integrate(fdx, bounds[0], bounds[1], _settings._angular_depth, 1.E-9, NULL)
                                             : gauss<double,N_GAUSS>::   integrate(fdx, bounds[0], bounds[1]);
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
                sum += kernel(previous->id(),s,t)*previous->basis_function(basis_id, t);
                sum *= this->_kinematics->jacobian(phi);
            };
            return sum;
        };
        return (_settings._adaptive_angular) ? gauss_kronrod<double,61>::integrate(fdx, _kinematics->phi_minus(s), _kinematics->phi_plus(s), _settings._angular_depth, 1.E-9, NULL)
                                             : gauss<double,N_GAUSS>::   integrate(fdx, _kinematics->phi_minus(s), _kinematics->phi_plus(s));
    };
}; // namespace iterateKT