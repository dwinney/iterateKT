// Output basis functions for the zero-th iteration (i.e. just omnes)
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "phase_shift.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "plotter.hpp"
#include "data_set.hpp"

void test_phases()
{
    using namespace iterateKT;

    phase_shift delta0("bern/phase_pipi_0.dat", 114., 1);
    phase_shift delta1("bern/phase_pipi_1.dat",  80., 1);
    phase_shift delta2("bern/phase_pipi_2.dat", 800., 0);

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    std::array<double,2> bounds = {4, 150};
    timer.start();

    plot p1 = plotter.new_plot();
    p1.set_legend(0.3, 0.7);
    p1.set_ranges(bounds, {-0.5, 2});
    p1.add_curve(bounds, [&](double s){ return delta0(s)/PI;}, solid(jpacColor::Blue,  "#it{S}0"));
    p1.add_curve(bounds, [&](double s){ return delta1(s)/PI;}, solid(jpacColor::Red,   "#it{P}1"));
    p1.add_curve(bounds, [&](double s){ return delta2(s)/PI;}, solid(jpacColor::Green, "#it{S}2"));
    p1.set_labels("#it{s} / m_{#pi}^{2}", "#delta(#it{s}) / #pi");
    p1.save("phases.pdf");

    timer.stop();
    timer.print_elapsed();
}; 