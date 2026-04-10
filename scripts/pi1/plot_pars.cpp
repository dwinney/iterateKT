// Take results from a multi-dimensional fit and output distributions of parameters
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
#include "data_set.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"
#include "COMPASS_pi1/fitter.hpp"
#include "COMPASS_pi1/data.hpp"

void plot_pars()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    // -----------------------------------------------------------------------
    // Operating options

    // Which range of m3pi bins to consider
    int min = 11, max = 49; 

    // Path to precalculated isoabrs
    std::string iso_path    = main_dir()+"/scripts/pi1/basis_functions/";
    // and the prefix given to each file
    std::string file_prefix = "CCD";

    // File containing parameters
    std::string path = main_dir()+"/scripts/pi1/pars/";

    // -----------------------------------------------------------------------
    // Import fit values

    auto minimal_alpha       = import_data<3>(path+"minimal_modNc.dat");
    auto minimal_moddelta    = import_data<3>(path+"minimal_modNd.dat");
    auto minimal_argdelta    = import_data<3>(path+"minimal_argNd.dat");
    
    auto nonminimal_alpha    = import_data<3>(path+"nonminimal_modNc.dat");
    auto nonminimal_moddelta = import_data<3>(path+"nonminimal_modNd.dat");
    auto nonminimal_argdelta = import_data<3>(path+"nonminimal_argNd.dat");

    auto nonminimal_beta     = import_data<3>(path+"nonminimal_modNcp.dat");

    // -----------------------------------------------------------------------
    // Plot results
    
    plotter plotter;

    // Plot distributions of parameters
    plot p1 = plotter.new_plot();
    p1.set_legend(0.675, 0.725);
    p1.add_horizontal(0);
    p1.set_labels("#it{m}_{3#pi}   [GeV]", "|#it{N}#kern[-0.5]{_{#it{c}}}| / 10^{3}");
    p1.add_data(nonminimal_alpha[0], {nonminimal_alpha[1]/1E3, nonminimal_alpha[2]/1E3},   dot(jpacColor::DarkGrey, "Non-minimal"));
    p1.add_data(   minimal_alpha[0],  {   minimal_alpha[1]/1E3,    minimal_alpha[2]/1E3},  dot(jpacColor::Blue,         "Minimal"));
    p1.save("Nc.pdf");

    plot p2 = plotter.new_plot();
    p2.set_legend(0.675, 0.725);
    p2.set_logscale(true);
    p2.set_labels("#it{m}_{3#pi}   [GeV]", "|#it{N}#kern[-0.5]{_{#it{d}}}| / 10^{3}");
    p2.add_data(nonminimal_moddelta[0], {nonminimal_moddelta[1]/1E3, nonminimal_moddelta[2]/1E3},  dot(jpacColor::DarkGrey,  "Non-minimal"));
    p2.add_data(   minimal_moddelta[0], {   minimal_moddelta[1]/1E3,    minimal_moddelta[2]/1E3},  dot(jpacColor::Red,           "Minimal"));
    p2.save("modNd.pdf");

    plot p3 = plotter.new_plot();
    p3.set_legend(0.675, 0.725);
    p3.add_horizontal(-PI, {2011, kDashed});
    p3.set_labels("#it{m}_{3#pi}   [GeV]", "#it{#phi}#kern[-0.1]{_{#it{d}}} - #it{#phi}#kern[-0.1]{_{#it{c}}}");
    p3.add_data(nonminimal_argdelta[0], {nonminimal_argdelta[1], nonminimal_argdelta[2]},  dot(jpacColor::DarkGrey,  "Non-minimal"));
    p3.add_data(   minimal_argdelta[0], {   minimal_argdelta[1],    minimal_argdelta[2]},  dot(jpacColor::Green,         "Minimal"));
    p3.save("argNd.pdf");

    plot p6 = plotter.new_plot();
    p6.set_legend(0.675, 0.725);
    p6.add_horizontal(0);
    p6.set_labels("#it{m}_{3#pi}   [GeV]", "|#it{N}#kern[-0.5]{_{#it{c}}}#kern[-0.6]{#it{'}}| / 10^{3}");
    p6.add_data(nonminimal_beta[0], {nonminimal_beta[1]/1E3, nonminimal_beta[2]/1E3},   dot(jpacColor::Orange));
    p6.save("Ncp.pdf");
};