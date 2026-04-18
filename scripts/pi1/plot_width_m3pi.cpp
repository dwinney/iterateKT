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

void plot_width_m3pi()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    // -----------------------------------------------------------------------
    // Operating options

    // Which range of m3pi bins to consider
    int min = 11, max = 49; 

    // Which tbin to plot
    std::array<uint,2> tbins = {0, 3};
    // IF to plot minimal or nonminimal
    bool minimal = true;
    
    // Path to precalculated isoabrs
    std::string iso_path    = main_dir()+"/analysis/COMPASS_pi1/basis_functions/";
    // and the prefix given to each file
    std::string file_prefix = "CCD";
    // path to par file
    std::string par_file    = (minimal) ? main_dir()+"/analysis/COMPASS_pi1/pars/minimal.dat"
                                        : main_dir()+"/analysis/COMPASS_pi1/pars/nonminimal.dat";

    // -----------------------------------------------------------------------
    // Data set up

    // Import our data sets  
    std::array<std::vector<data_set>,2> data;
    for (int j = 0; j < 2; j++)
    {
       for (int i = min; i <= max; i++) data[j].emplace_back(COMPASS::parse_JSON(i, tbins[j]));
    };

    // -----------------------------------------------------------------------
    // Import fit values

    std::vector<complex> minimal_pars    = COMPASS::import_parameters({min,max}, par_file);
    std::vector<complex> minimal_pars_c, minimal_pars_d;
    for (int i = 0; i < (minimal_pars.size()-3)/3; i++)
    {
        // Contact only
        minimal_pars_c.push_back( minimal_pars[3*i]   );
        minimal_pars_c.push_back( 0. );
        minimal_pars_c.push_back( 0. );
        // Deck only
        minimal_pars_d.push_back( 0. );
        minimal_pars_d.push_back( minimal_pars[3*i+1]  );
        minimal_pars_d.push_back( minimal_pars[3*i+2]   );
    };
    for (int n = 0; n < 3; n++)
    {
        int m = minimal_pars.size()-3;
        minimal_pars_c.push_back(minimal_pars[m+n]);
        minimal_pars_d.push_back(minimal_pars[m+n]);
    };

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    std::vector<double> m3pis = COMPASS::m3pi_bins();

    // Set up our amplitude (uniterated)
    amplitude amp  = new_amplitude<pi1_binned>(nullptr, std::make_tuple(m3pis, COMPASS::t_bins));
    amp->import_solution(iso_path+file_prefix);
    amp->precompute_dalitz(300);
    
    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    std::array<std::vector<double>,2> mws, ews, mws_c, mws_d;
    for (int j = 0; j < 2; j++)
    {
        for (auto bin : data[j])
        {
            double ew = 0, mw = 0, nw = 0;
            for (int i = 0; i < bin._z.size(); i++) ew += norm(bin._z[i])*bin._dx[i];
            ews[j].push_back(ew/1E4);
            
            COMPASS::fit_2D::process_parameters(minimal_pars, amp);
            amp->set_option(option::set_tbin,         bin._extras["t_bin"]);
            amp->set_option(option::set_mbin_COMPASS, bin._extras["m3pi_bin"]);
            mws[j].push_back(amp->width()/1E4);
            COMPASS::fit_2D::process_parameters(minimal_pars_c, amp);
            amp->set_option(option::set_tbin,         bin._extras["t_bin"]);
            amp->set_option(option::set_mbin_COMPASS, bin._extras["m3pi_bin"]);
            mws_c[j].push_back(amp->width()/1E4);
            COMPASS::fit_2D::process_parameters(minimal_pars_d, amp);
            amp->set_option(option::set_tbin,         bin._extras["t_bin"]);
            amp->set_option(option::set_mbin_COMPASS, bin._extras["m3pi_bin"]);
            mws_d[j].push_back(amp->width()/1E4);
        };
    };
        
        // -----------------------------------------------------------------------
        // Plot results

    plotter plotter;

    auto blue  = square(jpacColor::Blue,   "Full");
    blue._draw_opt = "L";
    auto red   = square(jpacColor::Red,    "Contact only");
    red._draw_opt = "PL";
    auto green = square(jpacColor::Green,  "Deck only");
    green._draw_opt = "PL";

    // Plot widths as a function of m3pi
    plot p1 = plotter.new_plot();
    p1.set_legend(0.65, 0.65);
    p1.add_header("#minus #it{t} = " + to_string(-COMPASS::t_bins[0], 2) + " GeV^{2}");
    p1.add_data(m3pis, mws_d[0],  green);
    p1.add_data(m3pis, mws_c[0],  red);
    p1.add_data(m3pis, mws[0],    blue); 
    p1.add_data(m3pis, ews[0],    dot(jpacColor::DarkGrey, "Data"));
    p1.set_labels("#it{m}_{3#pi}  [GeV]", "#Gamma(#it{t}, #it{m}_{3#pi}^{2}) / 10^{4}    [a.u.]");

    plot p2 = plotter.new_plot();
    p2.set_legend(0.65, 0.65);
    p2.add_header("#minus #it{t} = " + to_string(-COMPASS::t_bins[3], 2) + " GeV^{2}");
    p2.add_data(m3pis, mws_d[1],  green);
    p2.add_data(m3pis, mws_c[1],   red);
    p2.add_data(m3pis, mws[1],    blue);
    p2.add_data(m3pis, ews[1],    dot(jpacColor::DarkGrey, "Data"));
    p2.set_labels("#it{m}_{3#pi}  [GeV]", "#Gamma(#it{t}, #it{m}_{3#pi}^{2}) / 10^{4}    [a.u.]");

    plotter.combine({2,1}, {p1, p2}, "intensity_m3pi.pdf");
};