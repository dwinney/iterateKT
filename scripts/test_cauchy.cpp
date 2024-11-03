// Evaluate the functions which arise from cauchy integrals done analytically
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

void test_cauchy()
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
    settings._angular_integrator_depth    = 10;
    settings._dispersion_integrator_depth = 15;
    settings._interp_energy_low           = 2.;
    settings._interp_energy_high          = 20.;
    settings._interp_offset               = 0.2;
    settings._interp_points_low           = 100;
    settings._interp_points_high          = 50;
    settings._matching_interval           = 0.07;
    settings._expansion_eps               = 0.09;

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
    // p1.set_curve_points(200);
    p1.set_ranges({A,1.2}, {-2, 5});
    p1.print_to_terminal(true);
    p1.add_curve({A, 1.2}, [&](double s){ return std::real(first->regularized_integrand(0, s)); }, dashed(jpacColor::Blue));
    p1.add_curve({A, 1.2}, [&](double s){ return std::imag(first->regularized_integrand(0, s)); }, dashed(jpacColor::Red));
    p1.set_labels("#it{s} [GeV^{2}]", "disp");
    p1.add_vertical({A, C, D});
    p1.save("dispersion.pdf");

    timer.stop();
    timer.print_elapsed();
};