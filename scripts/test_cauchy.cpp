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
    using namespace V_to_3pi;

    // Set up general kinematics so everything knows masses
    kinematics kinematics = new_kinematics(M_OMEGA, M_PION);

    // Significant points in integration path
    double A = kinematics->A();
    double B = kinematics->B();
    double C = kinematics->C();
    double D = kinematics->D();

    // Set up our amplitude 
    amplitude amplitude = new_amplitude<isoscalar>(kinematics, "#Omega decay");

    // Just use the default settings
    // Define it outside so we can 
    settings settings = P_wave::default_settings();

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    amplitude->add_isobar<P_wave>(1, settings);
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
    iteration first = new_iteration(1, 3, grid, kinematics, settings);

    plot p1 = plotter.new_plot();
    p1.set_legend(0.45, 0.6);
    p1.set_curve_points(200);
    // p1.print_to_terminal(true);
    p1.add_curve({A, 1.}, [&](double s) { return std::real(first->regularized_integrand(0, s)); }, solid(jpacColor::Blue));
    p1.add_curve({A, 1.}, [&](double s) { return std::imag(first->regularized_integrand(0, s)); }, solid(jpacColor::Red));
    p1.set_labels("#it{s} [GeV^{2}]", "disp");
    p1.add_vertical({A, C, D});
    p1.add_horizontal(0);
    p1.save("dispersion.pdf");

    timer.stop();
    timer.print_elapsed();
};