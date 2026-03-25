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
    // Common central values for the 4 t-values using in all m3pi bins
    static const std::array<double,4>  t_bins = {-0.1205, -0.1675, -0.26, -0.663};

    // Also store the central m_bins for convenience, these are all 40 MeV in width
    static const std::array<double,39> m_bins = 
    {  0.96, 1.00, 1.04, 1.08, 1.12, 1.16, 1.20, 1.24, 1.28, 1.32, 1.36, 1.40, 1.44,
       1.48, 1.52, 1.56, 1.60, 1.64, 1.68, 1.72, 1.76, 1.80, 1.84, 1.88, 1.92, 1.96,
       2.00, 2.04, 2.08, 2.12, 2.16, 2.20, 2.24, 2.28, 2.32, 2.36, 2.40, 2.44, 2.48  };

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
            
        // Function to minimize is the average chi2 per dalitz plot
        static double fcn(std::vector<data_set> & data_vector, amplitude to_fit)
        {
            double chi2 = 0;
            for (auto data : data_vector)
            {
                to_fit->set_option(option::set_tbin,         data._extras["t_bin"]);
                to_fit->set_option(option::set_mbin_COMPASS, data._extras["m3pi_bin"]);
                for (int i = 0; i < data._N; i++)
                {
                    double from_data  = data._z[i];
                    double s = data._x[i], t = data._y[i], u = to_fit->get_kinematics()->Sigma() - s - t;
                    iterateKT::complex from_model = to_fit->evaluate(s, t, u);  

                    if (is_zero(data._dz[i])) continue;

                    // Minimize average chi2 per data plot
                    chi2  += norm((from_data - abs(from_model)) / data._dz[i]) / data._N; 
                };
            };
            return chi2 / data_vector.size();
        };
    };

    inline std::vector<complex> import_parameters(std::array<int,2> minmax, std::string in_pars_file)
    {
        std::vector<complex> pars;
        std::ifstream infile(in_pars_file);
        std::string line;

        int n = minmax[1] - minmax[0];

        if (!infile) warning("COMPASS::import_parameters", "Cannot open file " + in_pars_file + "!");

        int nimported = 0; // Mark how many lines we've imported
        while (std::getline(infile, line))
        {   
            if (line.empty())        continue; // skips empty lines
            if (line.front() == '#') continue; // Skip comment lines 
            std::istringstream is(line);  
            
            if (nimported < n+1)
            {
                double trash, alpha, redelta, imdelta;
                // Dont care about first two columns
                is >> trash >> trash;
                // we do about these though
                is >> alpha >> redelta >> imdelta;
    
                pars.push_back(alpha);
                pars.push_back(redelta+I*imdelta);
                nimported++;
                continue;
            };

            double b_alpha, b_delta;
            is >> b_alpha >> b_delta; 
            pars.push_back(b_alpha);
            pars.push_back(b_delta);
        };

        return pars;
    };

    inline void export_parameters(std::array<int,2> minmax, std::vector<complex> pars, std::string description, std::string out_pars_file)
    {    
        std::ofstream out;
        out.open(out_pars_file);
        int  precision = 12, spacing = precision + 10;

        int min = minmax[0], max = minmax[1];
        
        // Preamble info
        out << std::left << "# "+description << std::endl;
        out << std::left << "# "+std::string(5*spacing-2, '-') << std::endl;
        std::array<std::string,5> headers = {"# bin", "m3pi [GeV]", "alpha", "Re delta", "Im delta"};
        out << std::left;
        for (auto x : headers) out << std::setw(spacing) << x;
        out << endl;
        out << std::left << "# "+std::string(5*spacing-2, '-') << std::endl;
        // Table of subtraction pars
        for (uint i = 0; i <= max - min; i++)
        {
            out << std::left << std::setprecision(precision);
            out << std::setw(spacing) << i + min;
            out << std::setw(spacing) << m_bins[i+min-11];
            out << std::setw(spacing) << real(pars[2*i]);
            out << std::setw(spacing) << real(pars[2*i+1]);
            out << std::setw(spacing) << imag(pars[2*i+1]);
            out << std::endl;
        };
        // Tack on the two t-slopes at the end
        out << std::left << "# "+std::string(2*spacing-2, '-') << std::endl;
        out << std::left << std::setw(spacing) << "# b_alpha" << std::setw(spacing) << "b_delta" << std::endl;
        out << std::left << "# "+std::string(2*spacing-2, '-') << std::endl;
        out << std::left << std::setw(spacing) << real(pars[2*(max-min)+2]) << std::setw(spacing) << real(pars[2*(max-min)+3]) << std::endl;
        out.close();
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
                for (int j = 0; j < 2; j++) new_pars.push_back(pars[j]*csqrt(t_bins[i])*exp(pars[2+j]*(t_bins[i]-t_bins[0])));
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
                to_fit->set_option(option::set_mbin, i);
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