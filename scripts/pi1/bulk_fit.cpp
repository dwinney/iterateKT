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
    int min = 13, max = 49; 
    bool initial_from_file   = true;
    std::string in_pars_file = "/scripts/pi1/in_pars.dat";
    
    line();
    divider();
    std::cout << "OPTIONS:" << std::endl;
    line();
    ask_for_value<int>(min, 13, "Minimum m3π bin?");
    ask_for_value<int>(max, 49, "Maximum m3π bin?");
    ask_for_value<std::string>(in_pars_file, "/scripts/pi1/in_pars.dat", "File with initial values?");
    divider();
    line();

    // Fitter's stopping tolerance
    double tolerance = 0.001;

    // Path to precalculated isoabrs
    std::string iso_path    = main_dir()+"/analysis/COMPASS_pi1/basis_functions/";
    // and the prefix given to each file
    std::string file_prefix = "CD";

    // Are we taking initial values from file? if so which?
    
    // Where do we export the fit parameter values
    std::string out_pars_file  = main_dir()+"/scripts/pi1/out_pars.dat";
    // Put a file description at the beginning
    std::string description = "contact + deck";
    
    // -----------------------------------------------------------------------
    // Data set up

    // Import our data sets  
    std::vector<data_set>    data;          // Store the data
    std::vector<double>      m3pi_vals;     // m3pi values for each data_set
    std::vector<std::string> labels;        // parameter labels
    
    for (int i = min; i <= max; i++)
    {
        for (int j = 0; j < 4; j++) data.emplace_back(COMPASS::parse_JSON(i, j));
        m3pi_vals.push_back(data.back()._extras["m3pi"]);
        labels.push_back("alpha_"+to_string(i));
        labels.push_back("delta_"+to_string(i));
    };
    labels.push_back("b_alpha");
    labels.push_back("b_delta");

    // -----------------------------------------------------------------------
    // Import initial values
        
    std::vector<complex>   initial_vals;  // starting values for fitting
    if (initial_from_file) initial_vals = COMPASS::import_parameters({min, max}, main_dir()+in_pars_file);
    else                   initial_vals = std::vector<complex>(2*(max-min+1)+2, 1.);

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    // Set up our amplitude (uniterated)
    auto args     = std::make_tuple(m3pi_vals, COMPASS::t_bins);
    amplitude amp = new_amplitude<pi1_binned>(nullptr, args);
    amp->set_name("π₁ → 3π");

    // and import the pre-calculated isobars
    amp->import_solution(iso_path+file_prefix);

    // precalculate isobars in the decay region
    amp->precompute_dalitz(300);

    // -----------------------------------------------------------------------
    // Set up fitter

    fitter<amplitude,COMPASS::fit_2D> fitter(amp, "Combined");
    fitter.set_tolerance(tolerance/2*1E3);
    fitter.set_print_level(4);
    fitter.set_strategy(0);

    // Add all bins
    fitter.add_data(data);

    // Add three t-slopes in addition to three subtraction coeffs
    fitter.add_extra_parameters(2);

    fitter.set_parameter_labels(labels);
    // Fix alphas to all be real (and positive)
    for (int i = min; i <= max; i++)
    {
        fitter.fix_argument("alpha_"+to_string(i), 0.); 
    };
    // t-slopes as well
    fitter.make_real("b_alpha"); 
    fitter.make_real("b_delta"); 
    
    fitter.do_fit(initial_vals);

    // -----------------------------------------------------------------------
    // Print fit results to out_file
    COMPASS::export_parameters({min, max}, fitter.pars(), 
                               "average χ² per Dalitz: "+iterateKT::to_string(fitter.fcn()), 
                               out_pars_file);

    // -----------------------------------------------------------------------
    // Also plot summary of the chi2s
    
    std::array<std::vector<double>,4> chi2s;
    double chi2_dof;
    line();
    
    double total_chi2 = 0, total_N = 0, avg_chi2 = 0.;
    for (int j = 0; j < 4; j++)
    {
        for (int i = 0; i <= max-min; i++)
        {
            auto bin = data[4*i+j];
            double chi2 = 0;
            amp->set_option(option::set_tbin,         bin._extras["t_bin"]);
            amp->set_option(option::set_mbin_COMPASS, bin._extras["m3pi_bin"]);
            for (int i = 0; i < bin._N; i++)
            {
                double from_data  = bin._z[i];
                double s = bin._x[i], t = bin._y[i];
                complex from_model = amp->evaluate_in_dalitz(s, t);  
                
                if (is_zero(bin._dz[i])) continue;
                chi2  += norm((from_data - abs(from_model)) / bin._dz[i]); 
            };
            chi2s[j].push_back(chi2/bin._N);
            total_N    += bin._N;
            total_chi2 += chi2;
            avg_chi2   += chi2/bin._N;
        };
    };
    
    chi2_dof = total_chi2/(total_N - 3*(max-min+1)+2);
    avg_chi2 /= 4*chi2s[0].size();
    
    print<20>("True chi2/dof =",  chi2_dof);
    print<20>("Average chi2/N =", avg_chi2);
    line();

    // -----------------------------------------------------------------------
    // Plot results
    
    plotter plotter;

    double smin = m3pi_vals.front(), smax = m3pi_vals.back();

    // Plot distributions of chi2s
    plot p1 = plotter.new_plot();
    p1.set_labels("#it{m}_{3#pi}   [GeV]", "#chi^{2} / #it{n}_{#sigma}");
    p1.set_legend(0.65, 0.725);
    p1.add_horizontal(chi2_dof, {kBlack, kDashed});
    p1.add_data(m3pi_vals, chi2s[3], star(    jpacColor::Orange, "#minus #it{t} = 0.66 GeV^{2}"));
    p1.add_data(m3pi_vals, chi2s[2], triangle(jpacColor::Green,  "#minus #it{t} = 0.26 GeV^{2}"));
    p1.add_data(m3pi_vals, chi2s[1], square(  jpacColor::Red,    "#minus #it{t} = 0.17 GeV^{2}"));
    p1.add_data(m3pi_vals, chi2s[0], dot(     jpacColor::Blue,   "#minus #it{t} = 0.12 GeV^{2}"));
    p1.save("chi2s.pdf");
};