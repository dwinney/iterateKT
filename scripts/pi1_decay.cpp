// Twice subtracted KT amplitudes for omega decay with only P-wave in [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2006.01058
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "plotter.hpp"

#include "isobars/pi1.hpp"


void pi1_decay()
{
    using namespace iterateKT;

    // Set up general kinematics so everything knows masses
    // Use masses in units of pion mass
    kinematics kinematics = new_kinematics(1.60, M_PION);
    
    // Significant points in integration path
    double A = kinematics->A();
    double B = kinematics->B();
    double C = kinematics->C();
    double D = kinematics->D();

    // Set up our amplitude 
    amplitude amplitude = new_amplitude(kinematics);

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    amplitude->add_isobar<P_wave>(1, id::P_wave);
    isobar pwave = amplitude->get_isobar(id::P_wave);
    
    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    double smax = 3; uint N = 4;

    plot p1 = plotter.new_plot();
    p1.set_legend(0.2, 0.7);
    p1.set_curve_points(1000);
    p1.set_ranges({-1, smax}, {-4, 7.5});
    p1.set_labels("#it{s}   [GeV^{2}]", "F(#it{s} + #it{i}#epsilon)");
    p1.add_vertical({A, C});
    p1.add_horizontal(0);
    p1.add_curve( {-1, smax}, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, "#Omega_{1}");
    p1.add_dashed({-1, smax}, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });

    std::array<std::string,4> labels = {"1st", "2nd", "3rd", "4th"};
    for (int i = 1; i <= N; i++)
    {
        amplitude->iterate();
        p1.add_curve( {-15, smax}, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, labels[i-1]);
        p1.add_dashed({-15, smax}, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });
    };
    p1.save("pi1_isobars.pdf");

    timer.stop();
    timer.print_elapsed();
};