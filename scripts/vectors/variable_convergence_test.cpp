// Twice subtracted KT amplitudes for ω decay with only P-wave in [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
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

#include "isobars/vector.hpp"
#include "amplitudes/vector.hpp"

#include "plotter.hpp"

void variable_convergence_test()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // Set up general kinematics so everything knows masses
    // Use masses in units of GeV
    kinematics kinematics = new_kinematics(1.0, M_PION);
    
    // Significant points in integration path
    double sth = kinematics->sth(); // 4mpi^2
    double pth = kinematics->pth(); // (momega - mpi)^2
    double rth = kinematics->rth(); // (momega + mpi)^2

    // Set up our amplitude 
    amplitude amplitude = new_amplitude<vector_decay>(kinematics);

    int nmax = 6;
    settings sets = default_settings();
    sets._iteration_multiplier = 0.2;
    sets._iteration_slope      = (1-sets._iteration_multiplier)/nmax;

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    isobar pwave = amplitude->add_isobar<P_wave>({0, 1, 2}, 3, id::P_wave, "pwave", sets);
    
    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    std::array<double,2> bounds = {0., 1.};
    plot p1 = plotter.new_plot();
    p1.set_legend(0.25, 0.6);
    p1.set_curve_points(100);
    p1.set_ranges(bounds, {-5, 8.0});
    p1.set_labels("#it{s} [GeV^{2}]", "F_{a}(#it{s} + #it{i}#epsilon)");
    p1.add_vertical(rth);
    p1.shade_region({sth, pth});
    p1.add_horizontal(0);
    p1.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, "#Omega_{1}^{1}");
    p1.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });

    plot p2 = plotter.new_plot();
    p2.set_legend(0.25, 0.6);
    p2.set_curve_points(100);
    p2.set_ranges(bounds, {-2.5, 4.5});
    p2.set_labels("#it{s} [GeV^{2}]", "F_{b}(#it{s} + #it{i}#epsilon)");
    p2.add_vertical(rth);
    p2.shade_region({sth, pth});
    p2.add_horizontal(0);
    p2.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(1, s+IEPS)); },  "#it{s} #Omega_{1}^{1}");
    p2.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(1, s+IEPS)); });

    plot p3 = plotter.new_plot();
    p3.set_legend(0.25, 0.6);
    p3.set_curve_points(100);
    p3.set_ranges(bounds, {-2.5, 4.5});
    p3.set_labels("#it{s} [GeV^{2}]", "F_{c}(#it{s} + #it{i}#epsilon)");
    p3.add_vertical(rth);
    p3.shade_region({sth, pth});
    p3.add_horizontal(0);
    p3.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(2, s+IEPS)); },  "#it{s}^{2} #Omega_{1}^{1}");
    p3.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(2, s+IEPS)); });

    std::vector<std::string> labels = {"1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th"};
    for (int i = 1; i <= nmax + 4; i++)
    {
        amplitude->iterate();
        if (i >= nmax)
        {
            p1.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, "#it{n} = " + to_string(i));
            p1.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });
            p2.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(1, s+IEPS)); }, "#it{n} = " + to_string(i));
            p2.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(1, s+IEPS)); });
            p3.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(2, s+IEPS)); }, "#it{n} = " + to_string(i));
            p3.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(2, s+IEPS)); });
        };
        timer.lap("Iteration " + std::to_string(i));
    };

    plotter.combine({3,1}, {p1,p2,p3}, "convergence.pdf");

    
    timer.stop();
    timer.print_elapsed();
};