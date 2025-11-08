// Methods to interface iterateKT::fitter
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef COMPASS_FITTER_HPP
#define COMPASS_FITTER_HPP

#include "COMPASS_pi1/data.hpp"

namespace iterateKT { namespace COMPASS
{
    // common central values for the 4 t-values using in all m3pi bins
    static const std::array<double,4> t_bins = {-0.1205, -0.1675, -0.26, -0.663};

    // This fitter takes a single data set and fits to it
    struct fit_single_bin
    {
        static std::string data_type(int i)
        {
            if (i == kDalitz) return "Dalitz Plot";
            else return "ERROR!";
        };

        // Dont need any additional processing, just save
        static std::vector<complex> process_parameters(std::vector<complex> pars, amplitude to_fit)
        {
            to_fit->set_parameters(pars); 
            return pars; 
        };
            
        // Function to minimize
        // Filters whether we're looking at the real or imaginary parts 
        static double fcn(std::vector<data_set> & data_vector, amplitude to_fit)
        {
            double chi2 = 0;
            for (auto data : data_vector)
            {
                to_fit->set_option(option::set_tbin,            data._extras["t_bin"]);
                to_fit->set_option(option::set_m3pibin_COMPASS, data._extras["m3pi_bin"]);
                for (int i = 0; i < data._N; i++)
                {
                    double from_data  = data._z[i];
                    double s = data._x[i], t = data._y[i], u = to_fit->get_kinematics()->Sigma() - s - t;
                    iterateKT::complex from_model = to_fit->evaluate(s, t, u);  

                    if (is_zero(data._dz[i])) continue;
                    chi2  += norm((from_data - abs(from_model)) / data._dz[i]); 
                };
            };
            return chi2;
        };
    };

    // This one on the other hand assumes we are looking at multiple tbins and have couplings with explicit t-dependence
    struct fit_across_tbins
    {
        // None of these change
        static std::string data_type(int i)
        { return fit_single_bin::data_type(i); };
        static double fcn(std::vector<data_set> & data_vector, amplitude to_fit)
        { return fit_single_bin::fcn(data_vector, to_fit); };

        // Take in the parameters with the extra form factor slopes at the end
        // We assume pars.size() = 4
        // first two are the subtraction coefficients
        // last  two are the t-slopes
        static std::vector<complex> process_parameters(std::vector<complex> pars, amplitude to_fit)
        {
            for (int i = 0; i < 4; i++)
            {
                std::vector<complex> new_pars;
                // g -> g * t exp{b(t-t_0)}
                for (int j = 0; j < 2; j++) new_pars.push_back(pars[j]*t_bins[i]*exp(pars[2+j]*(t_bins[i]-t_bins[0])));
                to_fit->set_option(option::set_tbin, i);
                to_fit->set_parameters(new_pars);                
            };  
            return pars;
        };
    };

    // Here we fit in 2D (both binned in m3pi and t)
    struct fit_2D
    {
        // None of these change
        static std::string data_type(int i)
        { return fit_single_bin::data_type(i); };
        static double fcn(std::vector<data_set> & data_vector, amplitude to_fit)
        { return fit_single_bin::fcn(data_vector, to_fit); };

        static std::vector<complex> process_parameters(std::vector<complex> pars, amplitude to_fit)
        {
            // Calculate number of bins
            int N = (pars.size()-2)/2;
            complex b_cont = pars.end()[-2]; // Second to last par
            complex b_deck = pars.end()[-1]; // Last par

            // Cycle through m3pibins
            for (int i = 0; i < N; i++)
            {
                to_fit->set_option(option::set_m3pibin, i);
                // and through tbins
                for (int j = 0; j < 4; j++)
                {
                    to_fit->set_option(option::set_tbin, j);
                    complex cont = pars[2*i]  *t_bins[j]*exp(b_cont*(t_bins[j]-t_bins[0]));
                    complex deck = pars[2*i+1]*t_bins[j]*exp(b_deck*(t_bins[j]-t_bins[0]));
                    to_fit->set_parameters({cont, deck});
                };
            };
            return pars;
        };
    };
}; /* namespace COMPASS */ }; /* namespace iterateKT */

#endif /* COMPASS_FITTER_HPP */