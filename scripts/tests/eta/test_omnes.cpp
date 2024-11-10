// Output basis functions for the zero-th iteration (i.e. just omnes)
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
#include "decays/pseudoscalar.hpp"
#include "plotter.hpp"

void test_omnes()
{
    using namespace iterateKT;
    using namespace pseudoscalar;

    // Set up general kinematics so everything knows masses
    kinematics eta = new_kinematics(M_ETA/M_PION, 1.);

    // Set up our amplitude 
    amplitude amp_I1 = new_amplitude<I1_transition>(eta);
    amp_I1->add_isobar<I1_transition::S0_wave>(2);
    amp_I1->add_isobar<I1_transition::P1_wave>(1);
    amp_I1->add_isobar<I1_transition::S2_wave>(1);

    // Isolate our pwave
    isobar S0 = amp_I1->get_isobar(kI1_S0);
    isobar P1 = amp_I1->get_isobar(kI1_P1);
    isobar S2 = amp_I1->get_isobar(kI1_S2);

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    std::array<double,2> bounds = {-10, 100};
    timer.start();

    plot p1 = plotter.new_plot();
    p1.set_curve_points(4000);
    p1.set_legend(0.3, 0.8);
    p1.add_curve(bounds, [&](double s){ return std::real(S0->basis_function(0, s+IEPS));}, solid(jpacColor::Red));
    p1.add_curve(bounds, [&](double s){ return std::imag(S0->basis_function(0, s+IEPS));}, solid(jpacColor::Blue));
    p1.set_labels("#it{s} / m_{#pi}^{2}", "#Omega_{0}^{0}(#it{s} + #it{i}#epsilon)");

    timer.lap("plot 1");

    plot p2 = plotter.new_plot();
    p2.set_curve_points(4000);
    p2.set_legend(0.3, 0.8);
    p2.add_curve(bounds, [&](double s){ return std::real(S0->basis_function(1, s+IEPS));}, solid(jpacColor::Red));
    p2.add_curve(bounds, [&](double s){ return std::imag(S0->basis_function(1, s+IEPS));}, solid(jpacColor::Blue));
    p2.set_labels("#it{s} / m_{#pi}^{2}", "#it{s} #Omega_{0}^{0}(#it{s} + #it{i}#epsilon)");

    timer.lap("plot 2");

    plot p3 = plotter.new_plot();
    p3.set_legend(0.3, 0.8);
    p3.set_curve_points(4000);
    p3.add_curve(bounds, [&](double s){ return std::real(P1->basis_function(0, s+IEPS));}, solid(jpacColor::Red));
    p3.add_curve(bounds, [&](double s){ return std::imag(P1->basis_function(0, s+IEPS));}, solid(jpacColor::Blue));
    p3.set_labels("#it{s} / m_{#pi}^{2}", "#Omega_{1}^{1}(#it{s} + #it{i}#epsilon)");

    timer.lap("plot 3");

    plot p4 = plotter.new_plot();
    p4.set_curve_points(4000);
    p4.set_legend(0.3, 0.8);
    p4.add_curve(bounds, [&](double s){ return std::real(S2->basis_function(0, s+IEPS));}, solid(jpacColor::Red));
    p4.add_curve(bounds, [&](double s){ return std::imag(S2->basis_function(0, s+IEPS));}, solid(jpacColor::Blue));
    p4.set_labels("#it{s} / m_{#pi}^{2}", "#Omega_{0}^{2}(#it{s} + #it{i}#epsilon)");

    timer.lap("plot 2");

    // Save plots
    plotter.combine({2,2}, {p1,p2,p3,p4}, "omnes.pdf");

    timer.stop();
    timer.print_elapsed();
}; 