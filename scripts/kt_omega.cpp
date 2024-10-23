// Script to calculate the P-wave omega amplitude following [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// References:
// [1] - https://arxiv.org/abs/2006.01058
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "constants.hpp"
#include "decays/V_to_3pi.hpp"

#include "plotter.hpp"

using namespace iterateKT;

void kt_omega()
{
    using namespace iterateKT;

    // Set up general kinematics so everything knows masses
    kinematics omega = new_kinematics(M_OMEGA, M_PION);

    // Set up our amplitude 
    amplitude A(omega);

    // We need to load our amplitude with our isobars 
    A.add_isobar<V_to_3pi::P_wave>(1);

    // Test our pwave
    isobar pwave = A.get_isobar(V_to_3pi::kP_wave);

    double e = 0.1;
    // -----------------------------------------------------------------------
    
    plotter plotter;

    plot p = plotter.new_plot();
    p.add_curve({EPS, 1.0}, [&](double s){ return std::real(pwave->basis_function(0,s+IEPS));} );
    p.add_curve({EPS, 1.0}, [&](double s){ return std::imag(pwave->basis_function(0,s+IEPS));} );
    p.save("test.pdf");

    print("FIN");
}; 