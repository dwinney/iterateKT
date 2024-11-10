// Output basis functions for the zero-th iteration (i.e. just omnes)
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "phaseshifts/GKPY.hpp"
#include "plotter.hpp"

void test_phases()
{
    using namespace iterateKT;

    double s0 = M_PION*M_PION;

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();
    std::array<double,2> bounds = {4, 100};
    plot p1 = plotter.new_plot();
    p1.set_legend(0.3, 0.7);
    p1.set_ranges(bounds, {-0.5, 2});
    p1.add_curve(bounds, [&](double s){ return GKPY::phase_shift(0, 0, s0*s)/PI;}, solid(jpacColor::Blue,  "#it{S}0"));
    p1.add_curve(bounds, [&](double s){ return GKPY::phase_shift(1, 1, s0*s)/PI;}, solid(jpacColor::Red,   "#it{P}1"));
    p1.add_curve(bounds, [&](double s){ return GKPY::phase_shift(2, 0, s0*s)/PI;}, solid(jpacColor::Green, "#it{S}2"));
    p1.add_vertical(1.3*1.3/s0);
    p1.set_labels("#it{s} / m_{#pi}^{2}", "#delta(#it{s}) / #pi");

    p1.save("phases.pdf");

    timer.stop();
    timer.print_elapsed();
}; 