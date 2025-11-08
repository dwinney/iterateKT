// Isobars relevant for the decay of meson with JP = 1- into 3pi as in Ref. [1]
// Also allow isobars which are both C-odd and C-even. 
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] -  https://arxiv.org/abs/2006.01058
// ------------------------------------------------------------------------------

#ifndef PI_ONE_ISOBARS_HPP
#define PI_ONE_ISOBARS_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shift.hpp"

namespace iterateKT
{ 
    // Ids for all our isobars
    enum class id : unsigned int { P_wave, Contact, Deck };

    // The P-wave is the dominant isobar
    // In terms of individual isobars this is the only one we need
    class P_wave : public raw_isobar
    {
        // -----------------------------------------------------------------------
        public: 
        
        // Constructor 
        P_wave(isobar_args args) : raw_isobar(args){};

        inline uint    angular_momentum()   { return 1; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            if (iso_id != get_id()) return 0.; // Only interact with self
            complex  k  = _kinematics->kacser(s), kz = _kinematics->kz(s,t);
            return -3/2*(k*k - kz*kz); // We've multiplied by k^2 which is why singularity_power() = 2
        };
    };
}; // namespace iterateKT 

#endif // PI_ONE_ISOBARS_HPP