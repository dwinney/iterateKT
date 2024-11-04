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
    kinematics omega = new_kinematics(M_OMEGA, M_PION);

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
    p1.add_curve({EPS, 1.0}, [&](double s){ return std::real(pwave->basis_function(0,s+IEPS));}, solid(jpacColor::Blue, "+ieps"));
    p1.add_curve({EPS, 1.0}, [&](double s){ return std::imag(pwave->basis_function(0,s+IEPS));}, solid(jpacColor::Red));
    p1.add_curve({EPS, 1.0}, [&](double s){ return std::real(pwave->basis_function(0,s-IEPS));}, dashed(jpacColor::Blue,"-ieps"));
    p1.add_curve({EPS, 1.0}, [&](double s){ return std::imag(pwave->basis_function(0,s-IEPS));}, dashed(jpacColor::Red));
    p1.set_labels("#it{s} [GeV^{2}]", "F^{0}_{0} = #Omega");

    timer.lap("plot 1");

    // This prints the second basis function which will be s*Omega
    plot p2 = plotter.new_plot();
    p2.add_curve({EPS, 1.0}, [&](double s){ return std::real(pwave->basis_function(1,s+IEPS));}, solid(jpacColor::Blue, "+ieps"));
    p2.add_curve({EPS, 1.0}, [&](double s){ return std::imag(pwave->basis_function(1,s+IEPS));}, solid(jpacColor::Red));
    p2.add_curve({EPS, 1.0}, [&](double s){ return std::real(pwave->basis_function(1,s-IEPS));}, dashed(jpacColor::Blue,"-ieps"));
    p2.add_curve({EPS, 1.0}, [&](double s){ return std::imag(pwave->basis_function(1,s-IEPS));}, dashed(jpacColor::Red));
    p2.set_labels("#it{s} [GeV^{2}]", "F^{0}_{1} = #it{s} #Omega");

    timer.lap("plot 2");

    plot p3 = plotter.new_plot();
    p3.set_curve_points(200);
    p3.set_ranges({0,10}, {0,0.3});
    p3.set_legend(0.5,0.4);
    p3.add_curve({omega->sth(), 10.0}, [&](double s){ return pwave->LHC(s);});
    p3.set_labels("#it{s} [GeV^{2}]", "sin#delta(#it{s}) / |#Omega(#it{s})|");

    timer.lap("plot 3");

    // Save plots
    plotter.combine({3,1}, {p1,p2,p3}, "omnes.pdf");

    timer.stop();
    timer.print_elapsed();

    // This should produce an error and NaN
    print(pwave->basis_function(2, 0.5+IEPS));
}; 