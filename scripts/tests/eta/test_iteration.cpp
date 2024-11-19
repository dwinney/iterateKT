// Evaluate the pinocchio integration
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "basis.hpp"
#include "decays/pseudoscalar_dI1.hpp"

#include "plotter.hpp"


void test_iteration()
{
    using namespace iterateKT;
    using namespace pseudoscalar;

    // Set up general kinematics so everything knows masses
    // Assume masses are given in terms of pion mass
    kinematics kin = new_kinematics(M_ETA/M_PION, 1.);

    // Significant points in integration path
    double A = kin->A();
    double B = kin->B();
    double C = kin->C();
    double D = kin->D();

    // Set up our amplitude 
    amplitude amp_dI1 = new_amplitude<dI1_transition>(kin);
    amp_dI1->add_isobar<dI1_transition::S0_wave>(2); 
    amp_dI1->add_isobar<dI1_transition::P1_wave>(1); 
    amp_dI1->add_isobar<dI1_transition::S2_wave>(0);


    // -----------------------------------------------------------------------
    // s
    // Grab out isobars for plotting
    isobar S0 = amp_dI1->get_isobar(kdI1_S0);
    isobar P1 = amp_dI1->get_isobar(kdI1_P1);
    isobar S2 = amp_dI1->get_isobar(kdI1_S2);

    timer timer;
    timer.start();

    amp_dI1->iterate();
    timer.lap("first iteration");
    amp_dI1->iterate();
    timer.lap("second iteration");

    plotter plotter;
    double smin =  -10;
    double smax =   40;

    auto plot_basis = [&](isobar isobar, int i, std::string label)
    {
        plot p = plotter.new_plot();
        p.set_legend(false);
        p.set_curve_points(400);
        p.add_curve( {smin, smax}, [&](double s){ return std::real(isobar->basis_function(0, i, s+IEPS)); }, dashed(jpacColor::Blue));
        p.add_curve( {smin, smax}, [&](double s){ return std::imag(isobar->basis_function(0, i, s+IEPS)); }, dashed(jpacColor::Red));
        p.add_curve( {smin, smax}, [&](double s){ return std::real(isobar->basis_function(i, s+IEPS)); },    solid(jpacColor::Blue, "Real"));
        p.add_curve( {smin, smax}, [&](double s){ return std::imag(isobar->basis_function(i, s+IEPS)); },    solid(jpacColor::Red,  "Imaginary"));
        p.set_labels("#it{s} / #it{m}_{#pi}^{2}", label);
        return p;
    };

    plot f0a = plot_basis(S0, 0, "F_{0}^{#alpha}(s)");
    plot f1a = plot_basis(P1, 0, "F_{1}^{#alpha}(s)");
    plot f2a = plot_basis(S2, 0, "F_{2}^{#alpha}(s)");
    plot f0b = plot_basis(S0, 1, "F_{0}^{#beta}(s)");
    plot f1b = plot_basis(P1, 1, "F_{1}^{#beta}(s)");
    plot f2b = plot_basis(S2, 1, "F_{2}^{#beta}(s)");
    plot f0g = plot_basis(S0, 2, "F_{0}^{#gamma}(s)");
    plot f1g = plot_basis(P1, 2, "F_{1}^{#gamma}(s)");
    plot f2g = plot_basis(S2, 2, "F_{2}^{#gamma}(s)");

    plotter.combine({3,3}, {f0a, f1a, f2a, f0b, f1b, f2b, f0g, f1g, f2g}, "iterations.pdf");

    timer.stop();
    timer.print_elapsed();
};