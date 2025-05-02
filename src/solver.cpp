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
    void solver::export_solution(std::string prefix, uint precision)
    {
        uint spacing   = precision + 10;

        for (int i = 0; i < _isobars.size(); i++)
        {
            std::ofstream output;
            auto isobar = _isobars[i];
            std::string name = isobar->name();
            // As a precaution to unnamed isobars overriding the same file,  
            // append the index it appears with
            if (name == "isobar") name += to_string(i);
            output.open(prefix + "_" + name + ".dat");

            auto last_iter = isobar->get_iteration();
            uint N         = isobar->_subtractions->N_basis();
            for (auto s : isobar->_s_list)
            {
                output << std::setprecision(precision) << std::left << std::setw(spacing) << s;
                for (int j = 0; j < N; j++)
                {
                    complex ksf_disc = last_iter->ksf_inhomogeneity(j, s);
                    output << std::left << std::setw(spacing) << real(ksf_disc) 
                                        << std::setw(spacing) << imag(ksf_disc);
                }
                output << std::endl;
            }
            output.close();
        };
        return;
    };
};