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

void omega()
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

    // Just use the default settings
    // Define it outside so we can 
    settings settings = P_wave::default_settings();

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    amplitude->add_isobar<P_wave>(1);
    
    // Iterate once
    isobar pwave = amplitude->get_isobar(kP_wave);

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    plot p1 = plotter.new_plot();
    p1.set_legend(0.2, 0.7);
    p1.set_curve_points(100);
    p1.set_ranges({-10, 60}, {-3, 6.5});
    p1.set_labels("#it{s} / m_{#pi}^{2}", "F_{1}(#it{s} + #it{i}#epsilon)");
    p1.add_vertical({A, B, C, D});
    p1.add_horizontal(0);
    p1.add_curve( {-10, 60.}, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, "#Omega_{1}");
    p1.add_dashed({-10, 60.}, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });

    for (int i = 1; i < 5; i++)
    {
        amplitude->iterate();
        p1.add_curve( {-10, 60.}, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, std::to_string(i)+"th");
        p1.add_dashed({-10, 60.}, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });
    };

    p1.save("omega.pdf");

    timer.stop();
    timer.print_elapsed();
};