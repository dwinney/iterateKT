// Take results from a multi-dimensional fit and plot the integrated widths
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Instituto de Ciencias Nucleares (ICN)
//               Universidad Nacional Autonoma de Mexico (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "constants.hpp"
#include "plotter.hpp"
#include "fitter.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"
#include "COMPASS_pi1/fitter.hpp"
#include "COMPASS_pi1/data.hpp"

void plot_intensity_m3pi()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    // -----------------------------------------------------------------------
    // Operating options

    // Which range of m3pi bins to consider
    int min = 11, max = 49; 

    // Which tbin to plot
    uint tbin = 2;

    // If we have two terms or three
    bool minimal = false;

    // Path to precalculated isoabrs
    std::string iso_path    = main_dir()+"/scripts/pi1/basis_functions/";
    // and the prefix given to each file
    std::string file_prefix = "CCD";

    // File containing parameters
    std::string in_pars_file   = (minimal) ? main_dir()+"/scripts/pi1/minimal.dat" : main_dir()+"/scripts/pi1/nonminimal.dat";

    // -----------------------------------------------------------------------
    // Data set up

    // Import our data sets  
    std::vector<data_set> data;
    for (int i = min; i <= max; i++) data.emplace_back(COMPASS::parse_JSON(i, tbin));

    // -----------------------------------------------------------------------
    // Import fit values

    std::vector<complex> pars = COMPASS::import_parameters({min,max}, in_pars_file);

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    std::vector<double> m3pis = COMPASS::m3pi_bins();

    // Set up our amplitude (uniterated)
    amplitude amp  = new_amplitude<pi1_binned>(nullptr, std::make_tuple(m3pis, COMPASS::t_bins));
    amp->import_solution(iso_path+file_prefix);
    COMPASS::fit_2D::process_parameters(pars, amp);

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution

    std::vector<double> ws, ews;
    for (auto bin : data)
    {
        double ew = 0, bin_width = 0.04;
        for (auto M : bin._z)  ew += norm(M*bin_width);
        ews.push_back(ew);

        amp->set_option(option::set_tbin,         bin._extras["t_bin"]);
        amp->set_option(option::set_mbin_COMPASS, bin._extras["m3pi_bin"]);

        double w = amp->width();
        ws.push_back(w);
        print(bin._extras["m3pi_bin"], w, ew);
    };

    // -----------------------------------------------------------------------
    // Plot results

    plotter plotter;

    // Plot widths as a function of m3pi
    plot p1 = plotter.new_plot();
    p1.set_legend(0.65, 0.7);
    p1.add_header("#minus #it{t} = " + to_string(-COMPASS::t_bins[tbin], 2) + " GeV^{2}");
    p1.add_curve(m3pis,  ws,  solid(jpacColor::Blue, "Model"));
    p1.add_data (m3pis, ews,  dot(jpacColor::DarkGrey, "Data"));
    p1.set_labels("#it{m}_{3#pi}  [GeV]", "Integrated Intensity [a. u.]");
    p1.save("widths.pdf");
   
};