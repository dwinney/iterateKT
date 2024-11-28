// Top level class which defines the 1 -> 3 decay process. 
// The template feeds in all the relevant user info for a specific process.
// So far we require all equal mass particles
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "amplitude.hpp"

namespace iterateKT
{
    // -----------------------------------------------------------------------
    // Evaluating an amplitude just sums isobars and their associated prefactors
    // in each channel.
    complex raw_amplitude::evaluate(complex s, complex t, complex u)
    {
        complex result = 0;

        // S_CHANNEL
        for (auto f : _isobars)
        {
            complex term = prefactor_s(f->get_id(), s, t, u);
            if (is_zero(term)) continue;
            result += term * f->evaluate(s);
        };

        // T_CHANNEL
        for (auto f : _isobars)
        {
            complex term = prefactor_t(f->get_id(), s, t, u);
            if (is_zero(term)) continue;
            result += term * f->evaluate(t);
        };

        // U_CHANNEL
        for (auto f : _isobars)
        {
            complex term = prefactor_u(f->get_id(), s, t, u);
            if (is_zero(term)) continue;
            result += term * f->evaluate(u);
        };

        return result;
    };

    // -----------------------------------------------------------------------
    // Calculate one iteration of the KT equations for each of the isobars
    // which have been initialized
    void raw_amplitude::iterate()
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

}; // namespace iterateKT