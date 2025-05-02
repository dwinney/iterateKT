// Twice subtracted KT amplitudes for J/ψ decay with P- and F-waves in [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2304.09736
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "plotter.hpp"

#include "isobars/vector.hpp"
#include "amplitudes/vector.hpp"


void jpsi_decay()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // Set up general kinematics so everything knows masses
    // Use masses in units of GeV
    kinematics kinematics = new_kinematics(M_JPSI, M_PION);
    
    // Thresholds
    double sth = kinematics->sth(), pth = kinematics->pth(), rth = kinematics->rth();

    // Set up our amplitude 
    amplitude amplitude = new_amplitude<vector_decay>(kinematics);

    // Settings to adjust threshold from the defaults which are set for the omega
    settings jpsi_settings = default_settings();
    jpsi_settings._intermediate_energy  = 15.;
    jpsi_settings._cutoff               = 40.;
    jpsi_settings._interpolation_points = {800, 10, 300};

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    isobar pwave = amplitude->add_isobar<P_wave>(2, id::P_wave, "pwave", jpsi_settings);
    isobar fwave = amplitude->add_isobar<F_wave>(2, id::F_wave, "fwave", jpsi_settings);

    // Iterate 9 times
    amplitude->timed_iterate(9);
    
    // // -----------------------------------------------------------------------
    
    plotter plotter;

    std::array<double,2> bounds = {0., 3.};

    plot p1 = plotter.new_plot();
    p1.set_legend(0.25, 0.7);
    p1.set_curve_points(100);
    p1.set_ranges(bounds, {0, 180});
    p1.set_labels("#sqrt{#it{s}} [GeV]", "#delta_{3}(#it{s})  [#circ]");
    p1.add_curve( bounds, [&](double w) { return (180/PI)*fwave->phase_shift(w*w); }, solid(jpacColor::Green));

    plot p2 = plotter.new_plot();
    p2.set_legend(0.25, 0.7);
    p2.set_curve_points(1000);
    p2.set_ranges(bounds, {-5, 12});
    p2.set_labels("#sqrt{#it{s}} [GeV]", "#Omega_{3}(#it{s})");
    p2.add_curve( bounds, [&](double w) { return real(fwave->omnes(w*w+IEPS)); }, "Real");
    p2.add_curve( bounds, [&](double w) { return imag(fwave->omnes(w*w+IEPS)); }, "Imag");

    plotter.combine({2,1}, {p1,p2}, "fwave.pdf");
    
    auto plot_iterations = [&](int j, std::string label)
    {
        bounds = {0., 2.};
        // Output function to access real and imaginary parts of isobar
        auto rP = [&](int i, int j){return [j,pwave,i](double w){return real(pwave->basis_function(i, j, w*w+IEPS));}; };
        auto iP = [&](int i, int j){return [j,pwave,i](double w){return imag(pwave->basis_function(i, j, w*w+IEPS));}; };

        plot p1 = plotter.new_plot();
        p1.set_legend(0.8, 0.7);
        p1.set_curve_points(1000);
        p1.set_labels("#sqrt{#it{s}} [GeV]", "Re "+label+"(#it{s} + #it{i}#epsilon)");
        p1.add_curve( bounds, rP(1,j), "1^{st}");
        p1.add_curve( bounds, rP(3,j), "3^{rd}");
        p1.add_curve( bounds, rP(6,j), "6^{th}");
        p1.add_curve( bounds, rP(9,j), "9^{th}");

        plot p2 = plotter.new_plot();
        p2.set_legend(0.8, 0.7);
        p2.set_curve_points(1000);
        p2.set_labels("#sqrt{#it{s}} [GeV]", "Im "+label+"(#it{s} + #it{i}#epsilon)");
        p2.add_curve( bounds, iP(1,j), "1^{st}");
        p2.add_curve( bounds, iP(3,j), "3^{rd}");
        p2.add_curve( bounds, iP(6,j), "6^{th}");
        p2.add_curve( bounds, iP(9,j), "9^{th}");

        return std::array<plot,2>{p1,p2};
    };

    auto p_a = plot_iterations(0, "F_{a}"), p_b = plot_iterations(1, "F_{b}");
    plotter.combine({2,2}, {p_a[0], p_a[1], p_b[0], p_b[1]}, "iterations.pdf");
};