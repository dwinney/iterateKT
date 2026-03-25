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
    int min = 11, max = 36; 
    double tolerance = 0.01;

    // Path to precalculated isoabrs
    std::string iso_path    = main_dir()+"/scripts/pi1/basis_functions/";
    // and the prefix given to each file
    std::string file_prefix = "CD";

    // Are we taking initial values from file? if so which?
    bool initial_from_file   = true;
    std::string in_pars_file = main_dir()+"/scripts/pi1/in_fit.dat";
    
    // Where do we export the fit parameter values
    std::string out_pars_file  = main_dir()+"/scripts/pi1/out_pars.dat";
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

    if (initial_from_file) initial_vals = COMPASS::import_parameters({min, max}, in_pars_file);
    else                   initial_vals = std::vector<complex>(2*(max-min+1)+2, 1.);

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
    COMPASS::export_parameters({min, max}, fitter.pars(), 
                               "average χ² per Dalitz: "+iterateKT::to_string(fitter.fcn()), 
                               out_pars_file);
};