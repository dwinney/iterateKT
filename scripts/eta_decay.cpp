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
    solver.add_isobar<dI1_S0>({0,1}, 2, id::dI1_S0); 
    solver.add_isobar<dI1_P1>({0},   1, id::dI1_P1); 
    solver.add_isobar<dI1_S2>({},    1, id::dI1_S2);
    solver.add_isobar<dI0_P1>({0},   1, id::dI0_P1);
    solver.add_isobar<dI2_P1>({0},   1, id::dI2_P1); 
    solver.add_isobar<dI2_S2>({},    1, id::dI2_S2);

    // -----------------------------------------------------------------------
    // Iterate N times

    int N = 5;
    
    timer timer;
    timer.start();
    for (int i = 1; i <= N; i++)
    {
        solver.iterate();
        timer.lap("iteration " + std::to_string(i));
    }
    timer.stop();
    timer.print_elapsed();

    // -----------------------------------------------------------------------
    // Plot Results

    plotter plotter;
    double smin =  -10;
    double smax =  +75;

    std::array<std::string,7> labels = {"0th", "1st", "2nd", "3rd", "4th", "5th", "6th"};

    auto plot_basis = [&](isobar isobar, int i, std::string label)
    {
        plot p = plotter.new_plot();
        p.set_curve_points(1000);
        p.set_legend(false);
        p.set_labels("#it{s} / #it{m}_{#pi}^{2}", label);
       
        p.add_curve({smin, smax}, [&](double s){ return std::real(isobar->basis_function(1, i, s+IEPS)); }, dotted(jpacColor::Red));
        p.add_curve({smin, smax}, [&](double s){ return std::imag(isobar->basis_function(1, i, s+IEPS)); }, dotted(jpacColor::Blue));
        p.add_curve({smin, smax}, [&](double s){ return std::real(isobar->basis_function(2, i, s+IEPS)); }, dashed(jpacColor::Red));
        p.add_curve({smin, smax}, [&](double s){ return std::imag(isobar->basis_function(2, i, s+IEPS)); }, dashed(jpacColor::Blue));
        p.add_curve({smin, smax}, [&](double s){ return std::real(isobar->basis_function(5, i, s+IEPS)); }, solid(jpacColor::Red));
        p.add_curve({smin, smax}, [&](double s){ return std::imag(isobar->basis_function(5, i, s+IEPS)); }, solid(jpacColor::Blue));
        return p;
    };

    // Grab out isobars for plotting
    isobar dI1_S0 = solver.get_isobar(id::dI1_S0);
    isobar dI1_P1 = solver.get_isobar(id::dI1_P1);
    isobar dI1_S2 = solver.get_isobar(id::dI1_S2);
    isobar dI0_P1 = solver.get_isobar(id::dI0_P1);
    isobar dI2_P1 = solver.get_isobar(id::dI2_P1);
    isobar dI2_S2 = solver.get_isobar(id::dI2_S2);

    plot f0a = plot_basis(dI1_S0, 0, "F_{0}^{#alpha}(#it{s})");
    plot f1a = plot_basis(dI1_P1, 0, "F_{1}^{#alpha}(#it{s})");
    plot f2a = plot_basis(dI1_S2, 0, "F_{2}^{#alpha}(#it{s})");
    plot f0b = plot_basis(dI1_S0, 1, "F_{0}^{#beta}(#it{s})");
    plot f1b = plot_basis(dI1_P1, 1, "F_{1}^{#beta}(#it{s})");
    plot f2b = plot_basis(dI1_S2, 1, "F_{2}^{#beta}(#it{s})");
    plot f0g = plot_basis(dI1_S0, 2, "F_{0}^{#gamma}(#it{s})");
    plot f1g = plot_basis(dI1_P1, 2, "F_{1}^{#gamma}(#it{s})");
    plot f2g = plot_basis(dI1_S2, 2, "F_{2}^{#gamma}(#it{s})");
    plot g1e = plot_basis(dI0_P1, 3, "G_{1}^{#epsilon}(#it{s})");
    plot h1t = plot_basis(dI2_P1, 4, "H_{1}^{#theta}(#it{s})");
    plot h2t = plot_basis(dI2_S2, 4, "H_{2}^{#theta}(#it{s})");


    f0a.set_legend(0.2, 0.3);
    plotter.combine({3,4}, {f0a, f1a, f2a, f0b, f1b, f2b, f0g, f1g, f2g, g1e, h1t, h2t}, "eta_isobars.pdf");
};