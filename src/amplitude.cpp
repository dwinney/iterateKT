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

namespace iteratedOKT
{
    // -----------------------------------------------------------------------
    // Evaluating an amplitude just sums isobars and their associated prefactors
    // in each channel.
    complex amplitude::operator()(complex s, complex t, complex u)
    {
        complex result = 0;

        // S_CHANNEL
        for (auto f : _isobars)
        {
            complex term = f->prefactor(S_CHANNEL, s, t, u);
            if (is_zero(term)) continue;
            result += term * f->evaluate(s);
        };

        // T_CHANNEL
        for (auto f : _isobars)
        {
            complex term = f->prefactor(T_CHANNEL, s, t, u);
            if (is_zero(term)) continue;
            result += term * f->evaluate(t);
        };

        // U_CHANNEL
        for (auto f : _isobars)
        {
            complex term = f->prefactor(U_CHANNEL, s, t, u);
            if (is_zero(term)) continue;
            result += term * f->evaluate(u);
        };

        return result;
    };

}; // namespace iteratedOKT