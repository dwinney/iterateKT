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
#include "decays/pseudoscalar.hpp"

#include "plotter.hpp"


void test_F1_inhomogeneity()
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
    amplitude amp_I1 = new_amplitude<I1_transition>(kin);
    amp_I1->add_isobar<I1_transition::S0_wave>(2);
    amp_I1->add_isobar<I1_transition::P1_wave>(1);
    amp_I1->add_isobar<I1_transition::S2_wave>(0);

    std::vector<isobar> previous = amp_I1->get_isobars();
    
    // Isolate our pwave
    isobar S0 = amp_I1->get_isobar(kI1_S0);
    isobar P1 = amp_I1->get_isobar(kI1_P1);
    isobar S2 = amp_I1->get_isobar(kI1_S2);

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    auto grid = P1->calculate_next(previous);
    timer.lap("grid");

    // first nontrivial iteration
    iteration first = new_iteration(kin, grid, pseudoscalar::default_settings());

    // for (int i = 4; i <= 10; i++) print(i, first->integral(0, i+IEPS));
    // exit(1);

    double smax = 20;

    plot p1 = plotter.new_plot();
    p1.set_legend(0.25, 0.2);
    p1.set_curve_points(100);
    p1.add_curve( {A, smax}, [&](double s){ return std::real(first->half_regularized_integrand(0, s)); }, solid(jpacColor::Blue, "Real"));
    p1.add_curve( {A, smax}, [&](double s){ return std::imag(first->half_regularized_integrand(0, s)); }, solid(jpacColor::Red,  "Imaginary"));
    p1.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{1}^{#alpha,0} / #nu^{3}");

    plot p2 = plotter.new_plot();
    p2.set_legend(0.25, 0.2);
    p2.set_curve_points(100);
    p2.add_curve( {A, smax}, [&](double s){ return std::real(first->regularized_integrand(0, s)); }, solid(jpacColor::Green,  "Real"));
    p2.add_curve( {A, smax}, [&](double s){ return std::imag(first->regularized_integrand(0, s)); }, solid(jpacColor::Orange, "Imaginary"));
    p2.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{1}^{#alpha,0} / #kappa^{3}");

    plot p3 = plotter.new_plot();
    p3.set_legend(0.25, 0.2);
    p3.set_curve_points(400);
    p3.add_curve( {4.1, smax}, [&](double s){ return std::real(first->integral(1, s+IEPS)); }, solid(jpacColor::Green,  "Real"));
    p3.add_curve( {4.1, smax}, [&](double s){ return std::imag(first->integral(1, s+IEPS)); }, solid(jpacColor::Orange, "Imaginary"));
    p3.set_labels("#it{s} / #it{m}_{#pi}^{2}", "I_{1}^{#alpha,1}");

    plotter.combine({3,1}, {p1,p2,p3}, "F1_inhomogeneity.pdf");

    timer.stop();
    timer.print_elapsed();
};