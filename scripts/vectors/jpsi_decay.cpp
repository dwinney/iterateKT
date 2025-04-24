// Twice subtracted KT amplitudes for J/psi decay with P- and F-waves in [1]
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

#include "isobars/vector.hpp"
#include "amplitudes/vector.hpp"

#include "plotter.hpp"

void jpsi_decay()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // Set up general kinematics so everything knows masses
    // Use masses in units of GeV
    kinematics kinematics = new_kinematics(3.096, M_PION);
    
    // Thresholds
    double sth = kinematics->sth(), pth = kinematics->pth(), rth = kinematics->rth();

    // Set up our amplitude 
    amplitude amplitude = new_amplitude<vector_decay>(kinematics);

    // Settings to adjust threshold from the defaults which are set for the omega
    settings jpsi_settings = P_wave::default_settings();
    jpsi_settings._intermediate_energy  = 15.;
    jpsi_settings._cutoff               = 40.;
    jpsi_settings._interpolation_points = {800, 10, 300};


    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    amplitude->add_isobar<P_wave>(2, id::P_wave, "pwave", jpsi_settings);
    amplitude->add_isobar<F_wave>(2, id::F_wave, "fwave", jpsi_settings);

    isobar pwave = amplitude->get_isobar(id::P_wave), fwave = amplitude->get_isobar(id::F_wave);

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

    bounds = {0., 2.};
    plot p3 = plotter.new_plot();
    p3.set_legend(0.8, 0.7);
    p3.set_curve_points(1000);
    p3.set_ranges(bounds, {-3, 4.5});
    p3.set_labels("#sqrt{#it{s}} [GeV]", "Re #it{F}_{a}(#it{s})");
    p3.add_curve( bounds, [&](double w) { return real(pwave->basis_function(1, 0, w*w+IEPS)); }, "1^{st}");
    p3.add_curve( bounds, [&](double w) { return real(pwave->basis_function(3, 0, w*w+IEPS)); }, "3^{rd}");
    p3.add_curve( bounds, [&](double w) { return real(pwave->basis_function(6, 0, w*w+IEPS)); }, "6^{th}");
    p3.add_curve( bounds, [&](double w) { return real(pwave->basis_function(9, 0, w*w+IEPS)); }, "9^{th}");

    plot p4 = plotter.new_plot();
    p4.set_legend(0.8, 0.7);
    p4.set_curve_points(1000);
    p4.set_ranges(bounds, {0, 6.5});
    p4.set_labels("#sqrt{#it{s}} [GeV]", "Im #it{F}_{a}(#it{s})");
    p4.add_curve( bounds, [&](double w) { return imag(pwave->basis_function(1, 0, w*w+IEPS)); }, "1^{st}");
    p4.add_curve( bounds, [&](double w) { return imag(pwave->basis_function(3, 0, w*w+IEPS)); }, "3^{rd}");
    p4.add_curve( bounds, [&](double w) { return imag(pwave->basis_function(6, 0, w*w+IEPS)); }, "6^{th}");
    p4.add_curve( bounds, [&](double w) { return imag(pwave->basis_function(9, 0, w*w+IEPS)); }, "9^{th}");

    plot p5 = plotter.new_plot();
    p5.set_legend(0.8, 0.7);
    p5.set_curve_points(1000);
    p5.set_ranges(bounds, {-2, 4});
    p5.set_labels("#sqrt{#it{s}} [GeV]", "Re #it{F}_{b}(#it{s})");
    p5.add_curve( bounds, [&](double w) { return real(pwave->basis_function(1, 1, w*w+IEPS)); }, "1^{st}");
    p5.add_curve( bounds, [&](double w) { return real(pwave->basis_function(3, 1, w*w+IEPS)); }, "3^{rd}");
    p5.add_curve( bounds, [&](double w) { return real(pwave->basis_function(6, 1, w*w+IEPS)); }, "6^{th}");
    p5.add_curve( bounds, [&](double w) { return real(pwave->basis_function(9, 1, w*w+IEPS)); }, "9^{th}");

    plot p6 = plotter.new_plot();
    p6.set_legend(0.8, 0.7);
    p6.set_curve_points(1000);
    p6.set_ranges(bounds, {-1, 6});
    p6.set_labels("#sqrt{#it{s}} [GeV]", "Im #it{F}_{b}(#it{s})");
    p6.add_curve( bounds, [&](double w) { return imag(pwave->basis_function(1, 1, w*w+IEPS)); }, "1^{st}");
    p6.add_curve( bounds, [&](double w) { return imag(pwave->basis_function(3, 1, w*w+IEPS)); }, "3^{rd}");
    p6.add_curve( bounds, [&](double w) { return imag(pwave->basis_function(6, 1, w*w+IEPS)); }, "6^{th}");
    p6.add_curve( bounds, [&](double w) { return imag(pwave->basis_function(9, 1, w*w+IEPS)); }, "9^{th}");

    plotter.combine({2,2}, {p3,p4,p5,p6}, "iterations.pdf");
};