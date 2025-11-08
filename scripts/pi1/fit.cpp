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

void fit()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // -----------------------------------------------------------------------
    // Operating options

    // List of which m3pi bins to consider
    std::vector<int> m3pi_bins = { 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };

    // -----------------------------------------------------------------------
    // Data set up

    // Import our data sets  
    std::vector<data_set>    data;          // Store the data
    std::vector<double>      m3pi_vals;     // m3pi values for each data_set
    std::vector<std::string> labels;        // parameter labels
    std::vector<complex>     initial_vals;  // starting values for fitting

    for (int i = 0; i < m3pi_bins.size(); i++)
    {
        for (int j = 0; j < 4; j++) data.emplace_back(COMPASS::parse_JSON(m3pi_bins[i], j));
        m3pi_vals.push_back(data.back()._extras["m3pi"]);
        labels.push_back("alpha_"+to_string(i));
        labels.push_back("delta_"+to_string(i));
    };
    initial_vals.push_back(1985.7969);
    initial_vals.push_back(complex(-2400.1453161,-970.591590827));
    initial_vals.push_back(1985.7969);
    initial_vals.push_back(complex(-2400.1453161,-970.591590827));
    initial_vals.push_back(1985.7969);
    initial_vals.push_back(complex(-2400.1453161,-970.591590827));
    //
    initial_vals.push_back(1985.7969);
    initial_vals.push_back(complex(-2400.1453161,-970.591590827));
    initial_vals.push_back(1804.5645);
    initial_vals.push_back(complex(-2272.8218802,-1202.82385026));
    initial_vals.push_back(2209.3792);
    initial_vals.push_back(complex(-2684.78143239,-1099.14409278));
    initial_vals.push_back(2587.0256);
    initial_vals.push_back(complex(-3354.28079487,-1174.83527852));
    initial_vals.push_back(3220.4605);
    initial_vals.push_back(complex(-4228.71824002,-906.893350382));
    initial_vals.push_back(4009.3046);
    initial_vals.push_back(complex(-5044.64318525,-737.120567143));
    initial_vals.push_back(4393.5825);
    initial_vals.push_back(complex(-5755.63943158,-623.563984981));
    initial_vals.push_back(4525.749);
    initial_vals.push_back(complex(-5898.30651121,-625.965662493));
    initial_vals.push_back(4520.5132);
    initial_vals.push_back(complex(-5844.85341304,-669.359022922));

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    // Set up our amplitude 
    auto binning   = std::make_tuple(m3pi_vals, COMPASS::t_bins);
    amplitude amp  = new_amplitude<pi1_binned>(nullptr, binning);
    amp->set_name("π₁ → 3π");

    // -----------------------------------------------------------------------
    // Set up fitter

    fitter<COMPASS::fit_2D> fitter(amp, "Combined");
    fitter.set_tolerance(1E-5);
    fitter.set_print_level(3);

    // Add all bins
    fitter.add_data(data);

    // Add three t-slopes in addition to three subtraction coeffs
    fitter.add_extra_parameters(2);
    labels.push_back("b_alpha");
    labels.push_back("b_delta");
    initial_vals.push_back(8.59674732359);
    initial_vals.push_back(3.865448834);

    fitter.set_parameter_labels(labels);
    for (int i = 0; i < m3pi_bins.size(); i++) fitter.fix_argument("alpha_"+to_string(i), 0.); 
    fitter.make_real("b_alpha"); 
    fitter.make_real("b_delta"); 
    // run
    fitter.do_fit(initial_vals);
};