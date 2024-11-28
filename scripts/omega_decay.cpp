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
#include "decays/vector.hpp"

#include "plotter.hpp"

void omega_decay()
{
    using namespace iterateKT;
    using vector =  iterateKT::vector;

    // Set up general kinematics so everything knows masses
    // Use masses in units of pion mass
    kinematics kinematics = new_kinematics(M_OMEGA/M_PION, 1.);
    
    // Significant points in integration path
    double A = kinematics->A();
    double B = kinematics->B();
    double C = kinematics->C();
    double D = kinematics->D();

    // Set up our amplitude 
    amplitude amplitude = new_amplitude<vector>(kinematics);

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    amplitude->add_isobar<vector::P_wave>(2);
    isobar pwave = amplitude->get_isobar(id::P_wave);
    
    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    plot p1 = plotter.new_plot();
    p1.set_legend(0.2, 0.7);
    p1.set_curve_points(100);
    p1.set_ranges({-15, 70}, {-4, 6.5});
    p1.set_labels("#it{s} / m_{#pi}^{2}", "F_{a}(#it{s} + #it{i}#epsilon)");
    p1.add_vertical({A, C, D});
    p1.add_horizontal(0);
    p1.add_curve( {-15, 70.}, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, "#Omega_{1}");
    p1.add_dashed({-15, 70.}, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });

    plot p2 = plotter.new_plot();
    p2.set_legend(0.2, 0.7);
    p2.set_curve_points(100);
    p2.set_ranges({-15, 70}, {-125, 200});
    p2.set_labels("#it{s} / m_{#pi}^{2}", "F_{b}(#it{s} + #it{i}#epsilon)");
    p2.add_vertical({A, C, D});
    p2.add_horizontal(0);
    p2.add_curve( {-15, 70.}, [&](double s) { return std::real(pwave->basis_function(1, s+IEPS)); },  "#it{s} #Omega_{1}");
    p2.add_dashed({-15, 70.}, [&](double s) { return std::imag(pwave->basis_function(1, s+IEPS)); });

    std::array<std::string,4> labels = {"1st", "2nd", "3rd", "4th"};
    for (int i = 1; i <= 4; i++)
    {
        amplitude->iterate();
        p1.add_curve( {-15, 70.}, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, labels[i-1]);
        p1.add_dashed({-15, 70.}, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });
        p2.add_curve( {-15, 70.}, [&](double s) { return std::real(pwave->basis_function(1, s+IEPS)); }, labels[i-1]);
        p2.add_dashed({-15, 70.}, [&](double s) { return std::imag(pwave->basis_function(1, s+IEPS)); });
    };

    plotter.combine({2,1}, {p1,p2}, "omega_isobars.pdf");

    timer.stop();
    timer.print_elapsed();
};