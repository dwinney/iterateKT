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
#include "decays/vector.hpp"

#include "plotter.hpp"


void test_inhomogeneity()
{
    using namespace iterateKT;
    using vector =  iterateKT::vector;

    // Set up general kinematics so everything knows masses
    // Assume masses are given in terms of pion mass
    kinematics omega = new_kinematics(M_OMEGA/M_PION, 1.);

    // Significant points in integration path
    double A = omega->A();
    double B = omega->B();
    double C = omega->C();
    double D = omega->D();

    // Set up our amplitude 
    amplitude amplitude = new_amplitude<vector>(omega);

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    amplitude->add_isobar<vector::P_wave>(1);
    auto previous = amplitude->get_isobars();

    // Isolate our pwave
    isobar pwave = amplitude->get_isobar(kP_wave);

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    auto grid = pwave->calculate_next(previous);
    timer.lap("grid");

    // first nontrivial iteration
    iteration first = new_iteration(omega, grid, vector::P_wave::default_settings());

    double smax = 60;

    plot p1 = plotter.new_plot();
    p1.set_legend(0.45, 0.6);
    p1.set_curve_points(100);
    p1.set_ranges({A, smax}, {0, 0.015});
    p1.add_curve( {A, smax}, [&](double s){ return std::real(first->regularized_integrand(0, s)); }, solid(jpacColor::Blue, "Real"));
    p1.add_curve( {A, smax}, [&](double s){ return std::imag(first->regularized_integrand(0, s)); }, solid(jpacColor::Red,  "Imaginary"));
    p1.set_labels("#it{s}/#it{m}_{#pi}^{2}", "#tilde{F}^{1}_{0} / #kappa^{3}");
    p1.add_vertical({A, C, D});

    plot p3 = plotter.new_plot();
    p3.set_legend(0.45, 0.4);
    p3.set_curve_points(100);
    p3.set_ranges({A, smax}, {-2.0, 0.6});
    p3.add_curve( {A, smax}, [&](double s){ return std::real(first->half_regularized_integrand(0, s)); }, solid(jpacColor::Green,  "Real"));
    p3.add_curve( {A, smax}, [&](double s){ return std::imag(first->half_regularized_integrand(0, s)); }, solid(jpacColor::Orange, "Imaginary"));
    p3.set_labels("#it{s}/#it{m}_{#pi}^{2}", "#tilde{F}^{1}_{0} / #nu^{3}");
    p3.add_vertical({A, C, D});

    plotter.combine({2,1}, {p3,p1}, "inhomogeneity.pdf");

    timer.stop();
    timer.print_elapsed();
};