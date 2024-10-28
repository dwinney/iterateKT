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
#include "decays/V_to_3pi.hpp"

#include "plotter.hpp"

using namespace iterateKT;

void test_pinocchio()
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

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    amplitude.add_isobar<V_to_3pi::P_wave>(2);
    auto previous = amplitude.get_isobars();

    // Isolate our pwave
    isobar pwave = amplitude.get_isobar(kP_wave);

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    // Integral of the first badid function
    timer.start();
    plot p1 = plotter.new_plot();
    p1.set_curve_points(100);
    p1.set_ranges({A, 1.2}, {-0.2, 0.15});
    p1.set_legend(0.45, 0.3);
    p1.add_curve({A, 1.2}, [&](double s){ return std::real(pwave->pinocchio_integral(0, s, previous)); }, solid(jpacColor::Blue, "Real"));
    p1.add_curve({A, 1.2}, [&](double s){ return std::imag(pwave->pinocchio_integral(0, s, previous));},  solid(jpacColor::Red,  "Imaginary"));
    p1.set_labels("#it{s} [GeV^{2}]", "#kappa^{3} #LT(1-#it{z}^{2}) F_{0}^{0}(#it{s})#GT");
    p1.add_vertical({A, B, C, D});

    timer.lap("plot 1");

    plot p2 = plotter.new_plot();
    p2.set_curve_points(100);
    p2.set_ranges({A, 1.2}, {-0.01, 0.05});
    p2.set_legend(0.45, 0.6);
    p2.add_curve({A, 1.2}, [&](double s){ return std::real(pwave->pinocchio_integral(1, s, previous)); }, solid(jpacColor::Blue, "Real"));
    p2.add_curve({A, 1.2}, [&](double s){ return std::imag(pwave->pinocchio_integral(1, s, previous));},  solid(jpacColor::Red,  "Imaginary"));
    p2.set_labels("#it{s} [GeV^{2}]", "#kappa^{3} #LT(1-#it{z}^{2}) F_{1}^{0}(#it{s})#GT");
    p2.add_vertical({A, B, C, D});

    timer.lap("plot 2");

    plotter.combine({2,1}, {p1,p2}, "pinocchio.pdf");

    timer.stop();
    timer.print_elapsed();
};