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
    using V_to_3pi::kP_wave;

    // Set up general kinematics so everything knows masses
    kinematics omega = new_kinematics(M_OMEGA, M_PION);

    // Significant points in integration path
    double A = omega->A();
    double B = omega->B();
    double C = omega->C();
    double D = omega->D();


    // Set up our amplitude 
    amplitude amplitude(omega);

    settings settings;
    settings._angular_integrator_depth = 10;
    settings._interp_energy_low  = 2.;
    settings._interp_offset      = 0.2;
    settings._interp_points_low  = 100;
    settings._interp_points_high = 100;
    settings._matching_interval = 0.07;
    settings._expansion_eps  = 0.09;

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    amplitude.add_isobar<V_to_3pi::P_wave>(2, settings);
    auto previous = amplitude.get_isobars();

    // Isolate our pwave
    isobar pwave = amplitude.get_isobar(kP_wave);

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    auto grid = pwave->calculate_next(previous);
    timer.lap("grid");
    // first nontrivial iteration
    iteration first = new_iteration(2, 3, grid, omega, settings);

    plot p1 = plotter.new_plot();
    p1.set_legend(0.45, 0.6);
    p1.set_curve_points(200);
    p1.set_ranges({0,1.2}, {-0.1, 0.15});
    p1.add_curve({A, 1.2}, [&](double s){ return std::real(first->ksf_inhomogeneity(0, s)); }, solid(jpacColor::Blue, "Real"));
    p1.add_curve({A, 1.2}, [&](double s){ return std::imag(first->ksf_inhomogeneity(0, s)); }, solid(jpacColor::Red,  "Imaginary"));
    p1.set_labels("#it{s} [GeV^{2}]", "#tilde{F}_{0}");
    p1.add_vertical({A, D});

    plot p2 = plotter.new_plot();
    p2.set_legend(0.45, 0.6);
    p2.set_ranges({0,1.2}, {-0.01, 0.04});
    p2.add_curve({A, 1.2}, [&](double s){ return std::real(first->ksf_inhomogeneity(1, s)); }, solid(jpacColor::Blue, "Real"));
    p2.add_curve({A, 1.2}, [&](double s){ return std::imag(first->ksf_inhomogeneity(1, s)); }, solid(jpacColor::Red,  "Imaginary"));
    p2.set_labels("#it{s} [GeV^{2}]", "#tilde{F}_{1}");
    p2.add_vertical({A, D});

    plot p3 = plotter.new_plot();
    p3.color_offset(2);
    p3.set_legend(0.45, 0.6);
    p3.set_curve_points(200);
    p3.set_ranges({0,1.2}, {-0.3, 0.6});
    p3.add_curve({A, 1.2}, [&](double s){ return std::real(first->full_inhomogeneity(0, s)); }, solid(jpacColor::Blue, "Real"));
    p3.add_curve({A, 1.2}, [&](double s){ return std::imag(first->full_inhomogeneity(0, s)); }, solid(jpacColor::Red,  "Imaginary"));
    p3.set_labels("#it{s} [GeV^{2}]", "#tilde{F}_{0} / #nu^{3}");
    p3.add_vertical({A, D});

    plot p4 = plotter.new_plot();
    p4.color_offset(2);
    p4.set_legend(0.45, 0.6);
    p4.set_ranges({0,1.2}, {-0.01, 0.17});
    p4.add_curve({A, 1.2}, [&](double s){ return std::real(first->full_inhomogeneity(1, s)); }, solid(jpacColor::Blue, "Real"));
    p4.add_curve({A, 1.2}, [&](double s){ return std::imag(first->full_inhomogeneity(1, s)); }, solid(jpacColor::Red,  "Imaginary"));
    p4.set_labels("#it{s} [GeV^{2}]", "#tilde{F}_{1} / #nu^{3}");
    p4.add_vertical({A, D});

    plotter.combine({2,2}, {p1,p2,p3,p4}, "inhomogeneity.pdf");

    timer.stop();
    timer.print_elapsed();
};