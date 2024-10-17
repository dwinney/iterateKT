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

using namespace iteratedOKT;

void kt_omega()
{
    using namespace iteratedOKT;

    // Set up general kinematics so everything knows masses
    kinematics omega = new_kinematics(M_OMEGA, M_PION);

    // Set up our amplitude 
    amplitude A(omega);

    // We need to load our amplitude with our isobars 
    A.add_isobar<V_to_3pi::P_wave>();

    print("FIN");
}; 