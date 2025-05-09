// Basis functions for K -> 3pi following [1]
// These are calcualted in the isospin limit. The 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2403.17570
// ------------------------------------------------------------------------------

#include <filesystem>
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "basis.hpp"
#include "plotter.hpp"

#include "amplitudes/kaon.hpp"
#include "isobars/kaon.hpp"

void calculate_isobars()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // -----------------------------------------------------------------------
    
    // Isospin averaged Kaon and charged pion
    kinematics kin = new_kinematics(M_KAON_AVG, M_PION_PM);

    // Set up our amplitude 
    solver solver(kin);

    // Adjust expansion near pth
    settings sets = default_settings();
    sets._expansion_offsets   = {1E-3, 1E-2, 1E-2};

    // Add all the isobars, note the order they are added will be the order
    // the basis functions are generated
    std::vector<uint> empty = {};
    isobar F0 = solver.add_isobar<I1_S0>({0, 1, 2}, 2, id::dI1_I1_S0, "F0", sets);
    isobar F1 = solver.add_isobar<I1_P1>({1},       1, id::dI1_I1_P1, "F1", sets);
    isobar F2 = solver.add_isobar<I1_S2>(empty,     2, id::dI1_I1_S2, "F2", sets); 
    isobar G1 = solver.add_isobar<I0_P1>({1},       1, id::dI1_I0_P1, "G1", sets);
    isobar H1 = solver.add_isobar<I2_P1>({0, 1},    1, id::dI3_I2_P1, "H1", sets);
    isobar H2 = solver.add_isobar<I2_S2>(empty,     2, id::dI3_I2_S2, "H2", sets);

    // Iterate N times
    solver.timed_iterate(5);

    // -----------------------------------------------------------------------
    // Plot Results

    // Create a directory to place files
    std::string out_dir = main_dir()+"/scripts/kaon/basis_functions";
    std::filesystem::create_directory(out_dir);

    plotter plotter;
    double smin =  +0.0;
    double smax =  +0.5;

    auto plot_basis = [&](isobar isobar, int i, std::string func)
    {
        plot p = plotter.new_plot();
        p.set_curve_points(1000);
        p.set_legend(false);

        auto rF = [&](int j){return [j,isobar,i](double s){return real(isobar->basis_function(j, i, s+IEPS));}; };
        auto iF = [&](int j){return [j,isobar,i](double s){return imag(isobar->basis_function(j, i, s+IEPS));}; };

        std::string label = func+"(#it{s} + #it{i}#epsilon)";
        p.set_labels("#it{s}  [GeV^{2}]", label);
        p.shade_region({kin->sth(), kin->pth()});
        p.add_vertical(kin->rth());
        p.add_curve({smin, smax}, rF(1), dotted(jpacColor::Red));
        p.add_curve({smin, smax}, iF(1), dotted(jpacColor::Blue));
        p.add_curve({smin, smax}, rF(2), dashed(jpacColor::Red));
        p.add_curve({smin, smax}, iF(2), dashed(jpacColor::Blue));
        p.add_curve({smin, smax}, rF(5), solid(jpacColor::Red));
        p.add_curve({smin, smax}, iF(5), solid(jpacColor::Blue));
        return p;
    };

    std::vector<plot> F_plots;
    F_plots.emplace_back(plot_basis(F0, 0, "F_{0}^{#alpha}"));
    F_plots.emplace_back(plot_basis(F1, 0, "F_{1}^{#alpha}"));
    F_plots.emplace_back(plot_basis(F2, 0, "F_{2}^{#alpha}"));
    F_plots.emplace_back(plot_basis(F0, 1, "F_{0}^{#beta}"));
    F_plots.emplace_back(plot_basis(F1, 1, "F_{1}^{#beta}"));
    F_plots.emplace_back(plot_basis(F2, 1, "F_{2}^{#beta}"));
    F_plots.emplace_back(plot_basis(F0, 2, "F_{0}^{#gamma}"));
    F_plots.emplace_back(plot_basis(F1, 2, "F_{1}^{#gamma}"));
    F_plots.emplace_back(plot_basis(F2, 2, "F_{2}^{#gamma}"));
    F_plots.emplace_back(plot_basis(F0, 3, "F_{0}^{#zeta}"));
    F_plots.emplace_back(plot_basis(F1, 3, "F_{1}^{#zeta}"));
    F_plots.emplace_back(plot_basis(F2, 3, "F_{2}^{#zeta}"));
    plotter.combine({3,4}, F_plots, out_dir+"/Fs.pdf");

    plot pG = plot_basis(G1, 4, "G_{1}^{#eta}");
    pG.save(out_dir+"/G.pdf");

    std::vector<plot> H_plots;
    H_plots.emplace_back(plot_basis(H1, 5, "H_{1}^{#mu}"));
    H_plots.emplace_back(plot_basis(H2, 5, "H_{2}^{#mu}"));
    H_plots.emplace_back(plot_basis(H1, 6, "H_{1}^{#nu}"));
    H_plots.emplace_back(plot_basis(H2, 6, "H_{2}^{#nu}"));
    plotter.combine({2,2}, H_plots, out_dir+"/Hs.pdf");

    // Export the solution so it can be more easily recalled later
    solver.export_solution(out_dir+"/basis");
};