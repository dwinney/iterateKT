// Outpit P-wave basis functions for the zero-th iteration (i.e. just omnes)
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
#include "decays/V_to_3pi.hpp"

#include "plotter.hpp"

using namespace iterateKT;

void test_omnes()
{
    using namespace iterateKT;
    using namespace V_to_3pi;

    // Set up general kinematics so everything knows masses
    kinematics omega = new_kinematics(M_OMEGA/M_PION, 1.);

    // Set up our amplitude 
    amplitude A = new_amplitude<isoscalar>(omega);

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    A->add_isobar<P_wave>(2);

    // Isolate our pwave
    isobar pwave = A->get_isobar(kP_wave);

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    // first we print basis function 1 which should be just Omega
    plot p1 = plotter.new_plot();
    p1.set_legend(0.3, 0.8);
    p1.add_curve({EPS, 60}, [&](double s){ return std::real(pwave->omnes(s+IEPS));}, solid(jpacColor::Blue, "+ieps"));
    p1.add_curve({EPS, 60}, [&](double s){ return std::imag(pwave->omnes(s+IEPS));}, solid(jpacColor::Red));
    p1.add_curve({EPS, 60}, [&](double s){ return std::real(pwave->omnes(s-IEPS));}, dashed(jpacColor::Blue,"-ieps"));
    p1.add_curve({EPS, 60}, [&](double s){ return std::imag(pwave->omnes(s-IEPS));}, dashed(jpacColor::Red));
    p1.set_labels("#it{s} / m_{#pi}^{2}", "#Omega_{1}^{1}");

    timer.lap("plot 1");

    plot p3 = plotter.new_plot();
    p3.set_curve_points(200);
    p3.set_legend(0.5,0.4);
    p3.add_curve({omega->sth(), 250}, [&](double s){ return pwave->LHC(s);});
    p3.set_labels("#it{s} / m_{#pi}^{2}", "sin#delta(#it{s}) / |#Omega(#it{s})|");

    timer.lap("plot 2");

    // Save plots
    plotter.combine({2,1}, {p1,p3}, "omnes.pdf");

    timer.stop();
    timer.print_elapsed();
}; 