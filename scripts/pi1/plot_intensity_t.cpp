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

void plot_intensity_t()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    // -----------------------------------------------------------------------
    // Operating options

    // Which tbins to plot
    std::array<int,3> m3pibins = {12, 22, 36};

    // If we have two terms or three
    bool minimal = true;

    // Path to precalculated isoabrs
    std::string iso_path    = main_dir()+"/scripts/pi1/basis_functions/";
    // and the prefix given to each file
    std::string file_prefix = "CCD";

    // File containing parameters
    std::string in_pars_file   = (minimal) ? main_dir()+"/scripts/pi1/minimal.dat" : main_dir()+"/scripts/pi1/nonminimal.dat";

    // -----------------------------------------------------------------------
    // Data set up

    // Import our data sets  
    std::vector<double> ts = COMPASS::mt_bins(), twidths;
    std::array<std::vector<data_set>,3> data;
    for (int j = 0; j < 3; j++)
    {
        for (int i = 0; i < 4; i++)
        {
            data[j].emplace_back(COMPASS::parse_JSON(m3pibins[j], i));
            if (j==0) twidths.push_back(data[j].back()._extras["tbin_width"]);
        };
    };

    // -----------------------------------------------------------------------
    // Import fit values

    std::vector<complex> pars = COMPASS::import_parameters({11, 49}, in_pars_file);

    // Filter pars for indididual terms
    std::vector<complex> pars_c, pars_d;
    for (int n = 0; n < (pars.size()-3)/3; n++)
    {
        // Contact only
        pars_c.push_back(pars[3*n]);
        pars_c.push_back(0);
        pars_c.push_back(0);
        // Deck only
        pars_d.push_back(0);
        pars_d.push_back(0);
        pars_d.push_back(pars[3*n+2]);
    };
    for (int n = 0; n < 3; n++)
    {
        int m = pars.size()-3;
        pars_c.push_back(pars[m+n]);
        pars_d.push_back(pars[m+n]);
    };

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution

    // Set up our amplitude (uniterated)
    amplitude amp  = new_amplitude<pi1_binned>(nullptr, std::make_tuple(COMPASS::m3pi_bins(), COMPASS::t_bins));
    amp->import_solution(iso_path+file_prefix);
    COMPASS::fit_2D::process_parameters(pars, amp);

    // Contact only
    amplitude c_only = new_amplitude<pi1_binned>(nullptr, std::make_tuple(COMPASS::m3pi_bins(), COMPASS::t_bins));
    c_only->import_solution(iso_path+file_prefix);
    COMPASS::fit_2D::process_parameters(pars_c, c_only);

    // Deck only
    amplitude d_only  = new_amplitude<pi1_binned>(nullptr, std::make_tuple(COMPASS::m3pi_bins(), COMPASS::t_bins));
    d_only->import_solution(iso_path+file_prefix);
    COMPASS::fit_2D::process_parameters(pars_d, d_only);

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution

    std::array<std::vector<double>,3> ws, wcs, wds, ews;
    for (int j = 0; j < 3; j++)
    {
        for (auto bin : data[j])
        {
            double ew = 0, bin_width = 0.04;
            for (int i = 0; i < bin._z.size(); i++) ew += norm(bin._z[i])*bin._dx[i];
            ews[j].push_back(ew);
            
            amp->set_option(option::set_tbin,            bin._extras["t_bin"]);
            amp->set_option(option::set_mbin_COMPASS,    bin._extras["m3pi_bin"]);
            c_only->set_option(option::set_tbin,         bin._extras["t_bin"]);
            c_only->set_option(option::set_mbin_COMPASS, bin._extras["m3pi_bin"]);
            d_only->set_option(option::set_tbin,         bin._extras["t_bin"]);
            d_only->set_option(option::set_mbin_COMPASS, bin._extras["m3pi_bin"]);
            
            ws[j].push_back(amp->width());
            wcs[j].push_back(c_only->width());
            wds[j].push_back(d_only->width());
        };
    };

    // -----------------------------------------------------------------------
    // Plot results

    plotter plotter;

    // Plot widths as a function of m3pi
    plot p1 = plotter.new_plot();
    p1.set_legend(0.7, 0.1);
    p1.set_logscale(false, true);
    p1.add_header("#it{m}_{3#pi} = " + to_string(COMPASS::m_bins[m3pibins[0]-11]) + " GeV");
    p1.set_labels("#minus #it{t}  [GeV^{2}]", "Integrated Intensity  [a.u.]");
    p1.add_curve( ts,  ws[0],             solid(jpacColor::Blue,   "Full"));
    p1.add_curve( ts,  wcs[0],            solid(jpacColor::Red,    "Contact only"));
    p1.add_curve( ts,  wds[0],            solid(jpacColor::Green,  "Deck Only"));
    p1.add_data ({ts, twidths},  ews[0],  dot(jpacColor::DarkGrey, "Data"));

    plot p2 = plotter.new_plot();
    p2.set_legend(0.7, 0.075);
    p2.set_logscale(false, true);
    p2.add_header("#it{m}_{3#pi} = " + to_string(COMPASS::m_bins[m3pibins[1]-11]) + " GeV");
    p2.set_labels("#minus #it{t}  [GeV^{2}]", "Integrated Intensity  [a.u.]");
    p2.add_curve( ts,  ws[1],             solid(jpacColor::Blue));
    p2.add_curve( ts,  wcs[1],            solid(jpacColor::Red));
    p2.add_curve( ts,  wds[1],            solid(jpacColor::Green));
    p2.add_data ({ts, twidths},  ews[1],  dot(jpacColor::DarkGrey));

    plot p3 = plotter.new_plot();
    p3.set_legend(0.7, 0.20);
    p3.set_logscale(false, true);
    p3.add_header("#it{m}_{3#pi} = " + to_string(COMPASS::m_bins[m3pibins[2]-11]) + " GeV");
    p3.set_labels("#minus #it{t}  [GeV^{2}]", "Integrated Intensity  [a.u.]");
    p3.add_curve( ts,  ws[2],             solid(jpacColor::Blue));
    p3.add_curve( ts,  wcs[2],            solid(jpacColor::Red));
    p3.add_curve( ts,  wds[2],            solid(jpacColor::Green));
    p3.add_data ({ts, twidths},  ews[2],  dot(jpacColor::DarkGrey));

    plotter.stack({p1,p2,p3}, "widths.pdf");   
};