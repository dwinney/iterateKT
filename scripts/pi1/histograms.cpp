// Take results from a multi-dimensional fit and output distributions of chi2
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

void histograms()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    // -----------------------------------------------------------------------
    // Operating options

    // Which range of m3pi bins to consider
    int min = 15, max = 44; 

    // Path to precalculated isoabrs
    std::string iso_path    = "/scripts/pi1/basis_functions/";
    // and the prefix given to each file
    std::string file_prefix = "CD";

    // File containing parameters
    std::string in_pars_file   = "/scripts/pi1/in_pars.dat";

    // -----------------------------------------------------------------------
    // Data set up

    // Import our data sets  
    std::array<std::vector<data_set>,4> data;
    for (int j = 0; j < 4; j++)
    {
        for (int i = min; i <= max; i++) data[j].emplace_back(COMPASS::parse_JSON(i, j));
    };

    // -----------------------------------------------------------------------
    // Import fit values

    std::vector<complex> pars;

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

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    std::vector<double> m3pi_vals; 
    for (auto mbin : data[0]) m3pi_vals.push_back(mbin._extras["m3pi"]);

    // Set up our amplitude (uniterated)
    auto args   = std::make_tuple(m3pi_vals, COMPASS::t_bins);
    amplitude amp  = new_amplitude<pi1_binned>(nullptr, args);
    amp->set_name("π₁ → 3π");

    // Import the pre-calculated isobars
    amp->import_solution(iso_path+file_prefix);
    // and set parameters from above
    COMPASS::fit_2D::process_parameters(pars, amp);

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution

    std::array<std::array<std::vector<double>,2>,4> chi2s;
    double total_chi2 = 0, total_N = 0;
    for (int j = 0; j < 4; j++)
    {
        for (auto bin : data[j])
        {
            double chi2 = 0;
            amp->set_option(option::set_tbin,         bin._extras["t_bin"]);
            amp->set_option(option::set_mbin_COMPASS, bin._extras["m3pi_bin"]);
            for (int i = 0; i < bin._N; i++)
            {
                double from_data  = bin._z[i];
                double s = bin._x[i], t = bin._y[i], u = amp->get_kinematics()->Sigma() - s - t;
                complex from_model = amp->evaluate(s, t, u);  
                
                if (is_zero(bin._dz[i])) continue;
                chi2  += norm((from_data - abs(from_model)) / bin._dz[i]); 
            };
            chi2s[j][0].push_back(bin._extras["m3pi"]);
            chi2s[j][1].push_back(chi2/bin._N);
            total_N    += bin._N;
            total_chi2 += chi2;
        };
    };

    double chi2_dof = total_chi2/(total_N - pars.size());

    // -----------------------------------------------------------------------
    // Plot results

    plotter plotter;

    plot p = plotter.new_plot();
    p.set_labels("#it{m}_{3#pi}   [GeV]", "#chi^{2} / #it{N}");
    p.set_legend(0.6, 0.7);
    p.set_ranges({1.0, 2.4}, {1, 13});
    p.add_horizontal(chi2_dof, {kBlack, kDashed});
    p.add_data(chi2s[3][0], chi2s[3][1], star(    jpacColor::Orange, "#minus #it{t} = 0.66"));
    p.add_data(chi2s[2][0], chi2s[2][1], triangle(jpacColor::Green,  "#minus #it{t} = 0.26"));
    p.add_data(chi2s[1][0], chi2s[1][1], square(  jpacColor::Red,    "#minus #it{t} = 0.17"));
    p.add_data(chi2s[0][0], chi2s[0][1], dot(     jpacColor::Blue,   "#minus #it{t} = 0.12"));
    p.save("chi2s.pdf");

};