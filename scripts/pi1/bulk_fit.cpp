// Fit KT amplitudes for π1(1600) decay with only P-wave to data from [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://inspirehep.net/literature/1898933
// ------------------------------------------------------------------------------

#include <algorithm>
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "plotter.hpp"
#include "fitter.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"
#include "COMPASS_pi1/fitter.hpp"
#include "COMPASS_pi1/data.hpp"

void bulk_fit()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    // -----------------------------------------------------------------------
    // Operating options

    // Which range of m3pi bins to consider
    int min = 11, max = 49; 
    double tolerance = 0.0001;

    // Path to precalculated isoabrs
    std::string iso_path    = "/scripts/pi1/basis_functions/";
    // and the prefix given to each file
    std::string file_prefix = "CD";

    // Are we taking initial values from file? if so which?
    bool initial_from_file   = true;
    std::string in_pars_file = "/scripts/pi1/in_pars.dat";
    
    // Where do we export the fit parameter values
    std::string out_pars_file  = "/scripts/pi1/out_pars.dat";
    // Put a file description at the beginning
    std::string description = "deck + contact, no form factor";
    
    // -----------------------------------------------------------------------
    // Data set up

    // Import our data sets  
    std::vector<data_set>    data;          // Store the data
    std::vector<double>      m3pi_vals;     // m3pi values for each data_set
    std::vector<std::string> labels;        // parameter labels
    std::vector<complex>     initial_vals;  // starting values for fitting

    for (int i = min; i <= max; i++)
    {
        for (int j = 0; j < 4; j++) data.emplace_back(COMPASS::parse_JSON(i, j));
        m3pi_vals.push_back(data.back()._extras["m3pi"]);
        labels.push_back("alpha_"+to_string(i));
        labels.push_back("delta_"+to_string(i));
    };

    // -----------------------------------------------------------------------
    // Import initial values

    std::vector<complex> intial_vals;

    if (initial_from_file)
    {
        std::ifstream infile(main_dir()+in_pars_file);
        std::string line;

        int nimported = 0; // Mark how many lines we've imported
        while (std::getline(infile, line))
        {   
            if (line.empty())        continue; // skips empty lines
            if (line.front() == '#') continue; // Skip comment lines 
            std::istringstream is(line);  
            
            if (nimported < max-min+1)
            {
                double trash, alpha, redelta, imdelta;
                // Dont care about first two columns
                is >> trash >> trash;
                // we do about these though
                is >> alpha >> redelta >> imdelta;
    
                initial_vals.push_back(alpha);
                initial_vals.push_back(redelta+I*imdelta);
                nimported++;
                continue;
            };

            double b_alpha, b_delta;
            is >> b_alpha >> b_delta; 
            initial_vals.push_back(b_alpha);
            initial_vals.push_back(b_delta);
        };
    }
    else initial_vals = std::vector<complex>(2*(max-min+1)+2, 1.);

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    // Set up our amplitude (uniterated)
    auto args   = std::make_tuple(m3pi_vals, COMPASS::t_bins);
    amplitude amp  = new_amplitude<pi1_binned>(nullptr, args);
    amp->set_name("π₁ → 3π");

    // and import the pre-calculated isobars
    amp->import_solution(iso_path+file_prefix);

    // -----------------------------------------------------------------------
    // Set up fitter

    fitter<amplitude,COMPASS::fit_2D> fitter(amp, "Combined");
    fitter.set_tolerance(tolerance*1E3);
    fitter.set_print_level(4);
    fitter.set_strategy(0);

    // Add all bins
    fitter.add_data(data);

    // Add three t-slopes in addition to three subtraction coeffs
    fitter.add_extra_parameters(2);
    labels.push_back("b_alpha");
    labels.push_back("b_delta");

    fitter.set_parameter_labels(labels);
    // Fix alphas to all be real (and positive)
    for (int i = min; i <= max; i++) fitter.fix_argument("alpha_"+to_string(i), 0.); 
    // t-slopes as well
    fitter.make_real("b_alpha"); 
    fitter.make_real("b_delta"); 

    fitter.do_fit(initial_vals);
    
    // -----------------------------------------------------------------------
    // Print fit results to out_file
    
    auto pars = fitter.pars();
 
    std::ofstream out;
    out.open(main_dir()+out_pars_file);
    int  precision = 12, spacing = precision + 10;
    
    // Preamble info
    out << std::left << "# "+description << std::endl;
    out << std::left << "# average χ²/Dalitz = "+to_string(fitter.fcn()) << std::endl;
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
        out << std::setw(spacing) << m3pi_vals[i];
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