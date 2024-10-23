// Simple container class to pass around interpolations of basis amplitudes 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef BASIS_GRID_HPP
#define BASIS_GRID_HPP

#include "utilities.hpp"

namespace iterateKT
{
    struct basis_grid
    {
        basis_grid(){};

        std::vector<double> _s_list;
        std::vector<std::vector<double>> _re_list, _im_list;
    };
};

#endif //SETTINGS_HPP