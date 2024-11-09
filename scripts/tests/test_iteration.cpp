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

void test_iteration()
{
    using namespace iterateKT;
    using namespace V_to_3pi;

    // Set up general kinematics so everything knows masses
    kinematics kinematics = new_kinematics(M_OMEGA/M_PION, 1.);

    // Significant points in integration path
    double A = kinematics->A();
    double B = kinematics->B();
    double C = kinematics->C();
    double D = kinematics->D();

    // Set up our amplitude 
    amplitude amplitude = new_amplitude<isoscalar>(kinematics, "#Omega decay");

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    amplitude->add_isobar<P_wave>(2);

    // Isolate our pwave
    isobar pwave = amplitude->get_isobar(kP_wave);

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    plot p1 = plotter.new_plot();
    p1.set_legend(0.45, 0.6);
    p1.set_curve_points(1000);
    p1.set_ranges({0,60}, {-4, 7});
    p1.set_labels("#it{s} / m_{#pi}^{2}", "F_{0}");
    p1.add_vertical({A, C, D});
    p1.add_horizontal(0);
    
    plot p2 = plotter.new_plot();
    p2.set_legend(0.45, 0.6);
    p2.set_curve_points(1000);
    p2.set_ranges({0,60}, {-120, 200});
    p2.set_labels("#it{s} / m_{#pi}^{2}", "F_{1}");
    p2.add_vertical({A, C, D});
    p2.add_horizontal(0);

    p1.add_curve({0, 60}, [&](double s) { return std::real(pwave->basis_function(0, s+I*EPS)); }, dashed(jpacColor::Blue));
    p1.add_curve({0, 60}, [&](double s) { return std::imag(pwave->basis_function(0, s+I*EPS)); }, dashed(jpacColor::Red));
    p2.add_curve({0, 60}, [&](double s) { return std::real(pwave->basis_function(1, s+I*EPS)); }, dashed(jpacColor::Blue));
    p2.add_curve({0, 60}, [&](double s) { return std::imag(pwave->basis_function(1, s+I*EPS)); }, dashed(jpacColor::Red));
    amplitude->iterate();
    p1.add_curve({0, 60}, [&](double s) { return std::real(pwave->basis_function(0, s+I*EPS)); }, solid(jpacColor::Blue));
    p1.add_curve({0, 60}, [&](double s) { return std::imag(pwave->basis_function(0, s+I*EPS)); }, solid(jpacColor::Red));
    p2.add_curve({0, 60}, [&](double s) { return std::real(pwave->basis_function(1, s+I*EPS)); }, solid(jpacColor::Blue));
    p2.add_curve({0, 60}, [&](double s) { return std::imag(pwave->basis_function(1, s+I*EPS)); }, solid(jpacColor::Red));

    plotter.combine({2,1}, {p1,p2}, "dispersion.pdf");
    timer.stop();
    timer.print_elapsed();
};