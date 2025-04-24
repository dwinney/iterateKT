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
        for (int i = 0; i < _isobars.size(); i++) 
        {
           if  (_isobars[i]->calculate_inhomogeneity()) _isobars[i]->save_iteration(next[i]);
           else _isobars[i]->skip_iteration();
        };
        return;
    };

    // Iterate N times with nice little terminal messages
    void solver::timed_iterate(unsigned int N)
    {
        // Iterate
        divider();
        print("Solving KT with " + to_string(_subtractions->N_basis()) + " subtractions and " + to_string(N) + " iterations:");
        line();
        timer timer; 
        timer.start();
        for (int i = 0; i < N; i++)
        {
            iterate();
            timer.lap("iteration " + to_string(i+1));
        };
        timer.stop();
        line();

        timer.print_elapsed();
        divider();
    };

    // Print to file necessary info to reconstruct isobars later
    void solver::export_solution(std::string prefix)
    {
        for (int n = 0; n < _isobars.size(); n++)
        {
            std::ofstream output;
            
            auto isobar = _isobars[n];
            std::string name = isobar->name();
            // As a precaution to unnamed isobars overriding the same file,  
            // append the index it appears with
            if (name == "isobar") name += to_string(n);
            output.open(prefix + "_" + name + ".dat");

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