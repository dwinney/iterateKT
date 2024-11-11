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
#include "decays/vector.hpp"

#include "plotter.hpp"

void test_omnes()
{
    using namespace iterateKT;
    using vector =  iterateKT::vector;

    // Set up general kinematics so everything knows masses
    kinematics omega = new_kinematics(M_OMEGA/M_PION, 1.);

    // Set up our amplitude 
    amplitude A = new_amplitude<vector>(omega);

    // We need to load our amplitude with our isobars 
    // Up to two subtractions so we have two basis functions
    A->add_isobar<vector::P_wave>(2);

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
    p1.set_labels("#it{s} / m_{#pi}^{2}", "#Omega_{1}^{1}(#it{s})");

    p1.save("omnes.pdf");

    timer.stop();
    timer.print_elapsed();
}; 