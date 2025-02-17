// Minimal structure to solve the KT equations
// This should be used if only isobars are important and not any amplitude
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "solver.hpp"

namespace iterateKT
{
    // Calculate one iteration of the KT equations
    void solver::iterate()
    {
        if (_isobars.size() == 0)
        { warning("amplitude::iterate()", "No isobars have been initialized!"); return; }

        // To not advance to the next iteration before looping over all isobars
        //  we save the "next" discontinuities here first
        std::vector<basis_grid> next;

        // Each isobar takes full list of other isobars with which to calculate angular avgs
        for (auto previous : _isobars) next.emplace_back( previous->calculate_next(_isobars) );

        // Save all the new iterations thereby pushing every isobar up by one iteration
        for (int i = 0; i < next.size(); i++) _isobars[i]->save_iteration(next[i]);

        return;
    };

    // Print to file necessary info to reconstruct isobars later
    void solver::export_solution(std::string prefix)
    {
        for (auto isobar : _isobars)
        {
            std::ofstream output;
            output.open(prefix + "_" + isobar->name() + ".dat");

            auto grid = isobar->calculate_next(_isobars);
            for (int i = 0; i < grid._s_list.size(); i++)
            {
                output << std::left << std::setw(PRINT_SPACING) << grid._s_list[i];
                for (int j = 0; j < grid.N_basis(); j++)
                {
                    output << std::left << std::setw(PRINT_SPACING) << grid._re_list[j][i] 
                                        << std::setw(PRINT_SPACING) << grid._im_list[j][i];
                }
                output << std::endl;
            }
            output.close();
        };
        return;
    };
};