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
#include "basis_grid.hpp"
#include "decays/V_to_3pi.hpp"

#include "plotter.hpp"

using namespace iterateKT;

void test_inhomogeneity()
{
    using namespace iterateKT;
    using namespace V_to_3pi;

    // Set up general kinematics so everything knows masses
    kinematics omega = new_kinematics(M_OMEGA, M_PION);

    // Significant points in integration path
    double A = omega->A();
    double B = omega->B();
    double C = omega->C();
    double D = omega->D();


    // Set up our amplitude 
    amplitude amplitude = new_amplitude<isoscalar>(omega);

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    amplitude->add_isobar<V_to_3pi::P_wave>(2);
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
    iteration first = new_iteration(2, 3, grid, omega, P_wave::default_settings());

    plot p1 = plotter.new_plot();
    p1.set_legend(0.45, 0.6);
    p1.set_curve_points(200);
    p1.set_ranges({0,1.2}, {-0.1, 0.15});
    p1.add_curve({A, 1.2}, [&](double s){ return std::real(first->ksf_inhomogeneity(0, s)); }, solid(jpacColor::Blue, "Real"));
    p1.add_curve({A, 1.2}, [&](double s){ return std::imag(first->ksf_inhomogeneity(0, s)); }, solid(jpacColor::Red,  "Imaginary"));
    p1.set_labels("#it{s} [GeV^{2}]", "#tilde{F}^{1}_{0}");
    p1.add_vertical({A, C, D});

    plot p2 = plotter.new_plot();
    p2.set_legend(0.45, 0.6);
    p2.set_ranges({0,1.2}, {-0.01, 0.04});
    p2.add_curve({A, 1.2}, [&](double s){ return std::real(first->ksf_inhomogeneity(1, s)); }, solid(jpacColor::Blue, "Real"));
    p2.add_curve({A, 1.2}, [&](double s){ return std::imag(first->ksf_inhomogeneity(1, s)); }, solid(jpacColor::Red,  "Imaginary"));
    p2.set_labels("#it{s} [GeV^{2}]", "#tilde{F}^{1}_{1}");
    p2.add_vertical({A, C, D});

    plot p3 = plotter.new_plot();
    p3.set_legend(0.45, 0.6);
    p3.set_curve_points(200);
    p3.set_ranges({0,1.2}, {-0.3, 0.6});
    p3.add_curve({A, 1.2}, [&](double s){ return std::real(first->half_regularized_integrand(0, s)); }, solid(jpacColor::Green, "Real"));
    p3.add_curve({A, 1.2}, [&](double s){ return std::imag(first->half_regularized_integrand(0, s)); }, solid(jpacColor::Orange,  "Imaginary"));
    p3.set_labels("#it{s} [GeV^{2}]", "#tilde{F}^{1}_{0} / #nu^{3}");
    p3.add_vertical({A, C, D});

    plot p4 = plotter.new_plot();
    p4.set_legend(0.45, 0.6);
    p4.set_ranges({0,1.2}, {-0.01, 0.17});
    p4.add_curve({A, 1.2}, [&](double s){ return std::real(first->half_regularized_integrand(1, s)); }, solid(jpacColor::Green, "Real"));
    p4.add_curve({A, 1.2}, [&](double s){ return std::imag(first->half_regularized_integrand(1, s)); }, solid(jpacColor::Orange,  "Imaginary"));
    p4.set_labels("#it{s} [GeV^{2}]", "#tilde{F}^{1}_{1} / #nu^{3}");
    p4.add_vertical({A, C, D});

    plotter.combine({2,2}, {p1,p2,p3,p4}, "inhomogeneity.pdf");

    timer.stop();
    timer.print_elapsed();
};