// Fit KT amplitudes for π1(1600) decay with only P-wave in [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2212.11767
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "plotter.hpp"
#include "fitter.hpp"

#include "amplitudes/pi1.hpp"

void plot_deck()
{
    using namespace iterateKT;
    using iterateKT::complex;

    plotter plotter;
    std::array<double,2> bounds = {-0.5,2.5};
    
    plot p1 = plotter.new_plot();
    p1.set_curve_points(400);
    p1.set_labels("#sigma  [GeV^{2}]",  "H_{Triangle}(t, m_{3#pi}^{2} #; #sigma)");
    p1.add_curve(bounds, [&](double sig){ return real(deck::projection(-0.1, 1.1, sig)); }, solid( jpacColor::Blue,  "#it{m}_{3#pi} = 1.1 GeV"));
    p1.add_curve(bounds, [&](double sig){ return imag(deck::projection(-0.1, 1.1, sig)); }, dashed( jpacColor::Blue));
    p1.add_curve(bounds, [&](double sig){ return real(deck::projection(-0.1, 1.6, sig)); }, solid( jpacColor::Red,   "#it{m}_{3#pi} = 1.6 GeV"));
    p1.add_curve(bounds, [&](double sig){ return imag(deck::projection(-0.1, 1.6, sig)); }, dashed( jpacColor::Red));
    p1.add_curve(bounds, [&](double sig){ return real(deck::projection(-0.1, 2.0, sig)); }, solid( jpacColor::Green, "#it{m}_{3#pi} = 2.0 GeV"));
    p1.add_curve(bounds, [&](double sig){ return imag(deck::projection(-0.1, 2.0, sig)); }, dashed( jpacColor::Green));
    p1.set_ranges(bounds, {0., 3.5});
    p1.add_vertical(0);
    p1.add_header("#it{t} = #minus 0.1 GeV^{2}");
    p1.set_legend(0.6, 0.71);

    plot p2 = plotter.new_plot();
    p2.set_curve_points(400);
    p2.set_labels("#sigma  [GeV^{2}]",  "H_{Triangle}(t, m_{3#pi}^{2} #; #sigma)");
    p2.add_curve(bounds, [&](double sig){ return real(deck::projection(-0.1, 1.6, sig)); }, solid( jpacColor::Blue,  "#it{t} = #minus 0.1 GeV^{2}"));
    p2.add_curve(bounds, [&](double sig){ return imag(deck::projection(-0.1, 1.6, sig)); }, dashed( jpacColor::Blue));
    p2.add_curve(bounds, [&](double sig){ return real(deck::projection(-0.5, 1.6, sig)); }, solid( jpacColor::Red,   "#it{t} = #minus 0.5 GeV^{2}"));
    p2.add_curve(bounds, [&](double sig){ return imag(deck::projection(-0.5, 1.6, sig)); }, dashed( jpacColor::Red));
    p2.add_curve(bounds, [&](double sig){ return real(deck::projection(-1, 1.6, sig)); }, solid( jpacColor::Green, "#it{t} = #minus 1.0 GeV^{2}"));
    p2.add_curve(bounds, [&](double sig){ return imag(deck::projection(-1, 1.6, sig)); }, dashed( jpacColor::Green));
    p2.shade_region({norm(2*M_PION), norm(1.6-M_PION)});
    p2.set_ranges(bounds, {0., 3.5});
    p2.add_header("#it{m}_{3#pi} = 1.6 GeV");
    p2.set_legend(0.6, 0.71);

    plotter.combine({2,1}, {p1,p2}, "loop.pdf");
};