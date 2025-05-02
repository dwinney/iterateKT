// Basis functions for eta -> 3pi following [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2111.02417
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "basis.hpp"
#include "plotter.hpp"
#include "solver.hpp"

#include "isobars/eta.hpp"

void eta_decay()
{
    using namespace iterateKT;

    // -----------------------------------------------------------------------
    
    // Set up general kinematics so everything knows masses
    // Assume masses are given in terms of pion mass
    kinematics kin = new_kinematics(M_ETA/M_PION, 1.);

    // Set up our amplitude 
    solver solver(kin);
    
    // Add all the isobars, note the order they are added will be the order
    // the basis functions are generated
    std::vector<uint> empty = {}; // Pass empty to isobars with no sub polynomials
    isobar F0 = solver.add_isobar<dI1_S0>(2,        id::dI1_S0); 
    isobar F1 = solver.add_isobar<dI1_P1>(1,        id::dI1_P1); 
    isobar F2 = solver.add_isobar<dI1_S2>(empty, 1, id::dI1_S2);
    isobar G1 = solver.add_isobar<dI0_P1>(1,        id::dI0_P1);
    isobar H1 = solver.add_isobar<dI2_P1>(1,        id::dI2_P1); 
    isobar H2 = solver.add_isobar<dI2_S2>(empty, 1, id::dI2_S2);

    // Iterate N times
    int N = 5;
    solver.timed_iterate(N);

    // -----------------------------------------------------------------------
    // Plot Results

    plotter plotter;
    double smin =  -10;
    double smax =  +75;

    auto plot_basis = [&](isobar isobar, int i, std::string label)
    {
        plot p = plotter.new_plot();
        p.set_curve_points(1000);
        p.set_legend(false);
        p.shade_region({kin->sth(), kin->pth()});
        p.add_vertical(kin->rth());

        label += "(#it{s} + #it{i}#epsilon)";
        p.set_labels("#it{s} / #it{m}_{#pi}^{2}", label);

        auto rF = [&](int j){return [j,isobar,i](double s){return real(isobar->basis_function(j, i, s+IEPS));}; };
        auto iF = [&](int j){return [j,isobar,i](double s){return imag(isobar->basis_function(j, i, s+IEPS));}; };
       
        p.add_curve({smin, smax}, rF(1), dotted(jpacColor::Red));
        p.add_curve({smin, smax}, iF(1), dotted(jpacColor::Blue));
        p.add_curve({smin, smax}, rF(2), dashed(jpacColor::Red));
        p.add_curve({smin, smax}, iF(2), dashed(jpacColor::Blue));
        p.add_curve({smin, smax}, rF(5), solid(jpacColor::Red));
        p.add_curve({smin, smax}, iF(5), solid(jpacColor::Blue));
        return p;
    };

    std::vector<plot> plots;
    plots.emplace_back(plot_basis(F0, 0, "F_{0}^{#alpha}"));
    plots.emplace_back(plot_basis(F1, 0, "F_{1}^{#alpha}"));
    plots.emplace_back(plot_basis(F2, 0, "F_{2}^{#alpha}"));
    plots.emplace_back(plot_basis(F0, 1, "F_{0}^{#beta}"));
    plots.emplace_back(plot_basis(F1, 1, "F_{1}^{#beta}"));
    plots.emplace_back(plot_basis(F2, 1, "F_{2}^{#beta}"));
    plots.emplace_back(plot_basis(F0, 2, "F_{0}^{#gamma}"));
    plots.emplace_back(plot_basis(F1, 2, "F_{1}^{#gamma}"));
    plots.emplace_back(plot_basis(F2, 2, "F_{2}^{#gamma}"));
    plots.emplace_back(plot_basis(G1, 3, "G_{1}^{#epsilon}"));
    plots.emplace_back(plot_basis(H1, 4, "H_{1}^{#theta}"));
    plots.emplace_back(plot_basis(H2, 4, "H_{2}^{#theta}"));

    plotter.combine({3,4}, plots, "eta_isobars.pdf");
};