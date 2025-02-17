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

#include "isobars/omega.hpp"
#include "amplitudes/omega.hpp"

#include "plotter.hpp"

void omega_decay()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // Set up general kinematics so everything knows masses
    // Use masses in units of GeV
    kinematics kinematics = new_kinematics(M_OMEGA, M_PION);
    
    // Significant points in integration path
    double sth = kinematics->sth(); // 4mpi^2
    double pth = kinematics->pth(); // (momega - mpi)^2
    double rth = kinematics->rth(); // (momega + mpi)^2

    // Set up our amplitude 
    amplitude amplitude = new_amplitude<omega>(kinematics);

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    amplitude->add_isobar<P_wave>({0, 1}, 2, id::P_wave, "pwave");
    isobar pwave = amplitude->get_isobar(id::P_wave);

    // Subtraction coefficients
    complex a = 3.00E2;
    complex b = 2.88*exp(I*1.85);
    amplitude->set_parameters({a, b});
    
    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    std::array<double,2> bounds = {0., 1.};
    plot p1 = plotter.new_plot();
    p1.set_legend(0.25, 0.7);
    p1.set_curve_points(100);
    p1.set_ranges(bounds, {-4, 7.});
    p1.set_labels("#it{s} [GeV^{2}]", "F_{a}(#it{s} + #it{i}#epsilon)");
    p1.add_vertical(rth);
    p1.shade_region({sth, pth});
    p1.add_horizontal(0);
    p1.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, "#Omega_{1}^{1}");
    p1.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });

    plot p2 = plotter.new_plot();
    p2.set_legend(0.25, 0.7);
    p2.set_curve_points(100);
    p2.set_ranges(bounds, {-2.5, 4.5});
    p2.set_labels("#it{s} [GeV^{2}]", "F_{b}(#it{s} + #it{i}#epsilon)");
    p2.add_vertical(rth);
    p2.shade_region({sth, pth});
    p2.add_horizontal(0);
    p2.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(1, s+IEPS)); },  "#it{s} #Omega_{1}^{1}");
    p2.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(1, s+IEPS)); });

    std::array<std::string,4> labels = {"1st", "2nd", "3rd", "4th"};
    for (int i = 1; i <= 4; i++)
    {
        amplitude->iterate();
        p1.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, labels[i-1]);
        p1.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });
        p2.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(1, s+IEPS)); }, labels[i-1]);
        p2.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(1, s+IEPS)); });
        timer.lap("Iteration " + std::to_string(i));
    };

    plotter.combine({2,1}, {p1,p2}, "omega_isobars.pdf");

    line();
    print("Cal. width:", amplitude->width());
    print("Exp. width:", 8.49E-3*0.893);
    line();
    
    timer.stop();
    timer.print_elapsed();
};