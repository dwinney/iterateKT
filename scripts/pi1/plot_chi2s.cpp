// Take results from a multi-dimensional fit and output distributions of chi2/N 
// for each dalitz plot as well as the fit parameters as a function of m3π
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

void plot_chi2s()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    // -----------------------------------------------------------------------
    // Operating options

    // Which range of m3pi bins to consider
    int min = 11, max = 49; 

    // Path to precalculated isoabrs
    std::string iso_path    = main_dir()+"/analysis/COMPASS_pi1/basis_functions/";
    // and the prefix given to each file
    std::string file_prefix = "CCD";

    // File containing parameters
    std::string file_path   = main_dir()+"/analysis/COMPASS_pi1/pars/";

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
    
    std::array<std::vector<complex>,2> pars;
    std::array<std::string,2> files = {file_path+"minimal.dat", file_path+"nonminimal.dat"};
    for (int i = 0; i < 2; i++) pars[i] = COMPASS::import_parameters({min,max}, files[i]);
    
   
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
    amp->precompute_dalitz(300);
    
    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    std::array<std::string,2> header = {"Minimal fit:", "Non-minimal fit:"};
    std::array<std::array<std::vector<double>,4>,2> chi2s;
    std::array<double,2> chi2_dof;
    line();
    for (int k = 0; k < 2; k++)
    {
        COMPASS::fit_2D::process_parameters(pars[k], amp);
        
        double total_chi2 = 0, total_N = 0, avg_chi2 = 0.;
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
                    double s = bin._x[i], t = bin._y[i];
                    complex from_model = amp->evaluate_in_dalitz(s, t);  
                    
                    if (is_zero(bin._dz[i])) continue;
                    chi2  += norm((from_data - abs(from_model)) / bin._dz[i]); 
                };
                chi2s[k][j].push_back(chi2/bin._N);
                total_N    += bin._N;
                total_chi2 += chi2;
                avg_chi2   += chi2/bin._N;
            };
        };
        
        chi2_dof[k] = total_chi2/(total_N - (3+k)*(max-min+1)+(2+k));
        avg_chi2 /= 4*chi2s[k][0].size();
        
        print<20>(header[k]);
        print<20>("True chi2/dof =",  chi2_dof[k]);
        print<20>("Average chi2/N =", avg_chi2);
        line();
    };

    // -----------------------------------------------------------------------
    // Plot results
    
    plotter plotter;

    double smin = m3pi_vals.front(), smax = m3pi_vals.back();

    // Plot distributions of chi2s
    plot p1 = plotter.new_plot();
    p1.set_labels("#it{m}_{3#pi}   [GeV]", "#chi^{2} / #it{n}_{#sigma}");
    p1.set_legend(0.65, 0.725);
    p1.set_ranges({0.9,2.6}, {0, 22});
    p1.add_horizontal(chi2_dof[0], {kBlack, kDashed});
    p1.add_data(m3pi_vals, chi2s[0][3], star(    jpacColor::Orange, "#minus #it{t} = 0.66 GeV^{2}"));
    p1.add_data(m3pi_vals, chi2s[0][2], triangle(jpacColor::Green,  "#minus #it{t} = 0.26 GeV^{2}"));
    p1.add_data(m3pi_vals, chi2s[0][1], square(  jpacColor::Red,    "#minus #it{t} = 0.17 GeV^{2}"));
    p1.add_data(m3pi_vals, chi2s[0][0], dot(     jpacColor::Blue,   "#minus #it{t} = 0.12 GeV^{2}"));
    p1.save("minimal_chi2s.pdf");

    // Plot distributions of chi2s
    plot p2 = plotter.new_plot();
    p2.set_labels("#it{m}_{3#pi}   [GeV]", "#chi^{2} / #it{n}_{#sigma}");
    p2.set_legend(0.65, 0.725);
    p2.set_ranges({0.9,2.6}, {0, 22});
    p2.add_horizontal(chi2_dof[1], {kBlack, kDashed});
    p2.add_data(m3pi_vals, chi2s[1][3], star(    jpacColor::Orange, "#minus #it{t} = 0.66 GeV^{2}"));
    p2.add_data(m3pi_vals, chi2s[1][2], triangle(jpacColor::Green,  "#minus #it{t} = 0.26 GeV^{2}"));
    p2.add_data(m3pi_vals, chi2s[1][1], square(  jpacColor::Red,    "#minus #it{t} = 0.17 GeV^{2}"));
    p2.add_data(m3pi_vals, chi2s[1][0], dot(     jpacColor::Blue,   "#minus #it{t} = 0.12 GeV^{2}"));
    p2.save("nonminimal_chi2s.pdf");
};