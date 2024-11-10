// Evaluate the pinocchio integration for C-conserving eta -> 3pi
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
#include "basis_grid.hpp"
#include "decays/pseudoscalar.hpp"

#include "plotter.hpp"


void test_inhomogeneity()
{
    using namespace iterateKT;
    using namespace pseudoscalar;

    // Set up general kinematics so everything knows masses
    // Assume masses are given in terms of pion mass
    kinematics eta = new_kinematics(M_ETA/M_PION, 1.);

    // Significant points in integration path
    double A = eta->A();
    double B = eta->B();
    double C = eta->C();
    double D = eta->D();

    // Set up our amplitude 
    amplitude amp_I1 = new_amplitude<I1_transition>(eta);
    amp_I1->add_isobar<I1_transition::S0_wave>(2);
    amp_I1->add_isobar<I1_transition::P1_wave>(1);
    amp_I1->add_isobar<I1_transition::S2_wave>(1);

    std::vector<isobar> previous = amp_I1->get_isobar();
    
    // Isolate our pwave
    isobar pwave = amp_I1->get_isobar(kI1_P1);

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    auto grid = pwave->calculate_next(previous);
    timer.lap("grid");

    // first nontrivial iteration
    iteration first = new_iteration(1, 3, grid, omega, vector::P_wave::default_settings());

    double smax = 60;

    plot p1 = plotter.new_plot();
    p1.set_legend(0.45, 0.6);
    p1.set_curve_points(100);
    p1.set_ranges({A, smax}, {-100, 100});
    p1.add_curve( {A, smax}, [&](double s){ return std::real(first->ksf_inhomogeneity(0, s)); }, solid(jpacColor::Blue, "Real"));
    p1.add_curve( {A, smax}, [&](double s){ return std::imag(first->ksf_inhomogeneity(0, s)); }, solid(jpacColor::Red,  "Imaginary"));
    p1.set_labels("#it{s}/#it{m}_{#pi}^{2}", "#tilde{F}^{1}_{0}");
    p1.add_vertical({A, C, D});

    plot p3 = plotter.new_plot();
    p3.set_legend(0.45, 0.4);
    p3.set_curve_points(100);
    p3.set_ranges({A, smax}, {-2.0, 0.6});
    p3.add_curve( {A, smax}, [&](double s){ return std::real(first->half_regularized_integrand(0, s)); }, solid(jpacColor::Green,  "Real"));
    p3.add_curve( {A, smax}, [&](double s){ return std::imag(first->half_regularized_integrand(0, s)); }, solid(jpacColor::Orange, "Imaginary"));
    p3.set_labels("#it{s}/#it{m}_{#pi}^{2}", "#tilde{F}^{1}_{0} / #nu^{3}");
    p3.add_vertical({A, C, D});

    plotter.combine({2,1}, {p1,p3}, "inhomogeneity.pdf");

    timer.stop();
    timer.print_elapsed();
};