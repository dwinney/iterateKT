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
#include "GKPY.hpp"
#include "phase_shift.hpp"

namespace iterateKT
{ 
    // Ids for all our isobars, we only have one though
    enum class id : unsigned int { P_wave };

    inline settings default_settings()
    {
        settings sets;
        sets._exclusion_points        = 6;
        sets._infinitesimal           = 1E-7;
        sets._intermediate_energy     = 4;
        sets._cutoff                  = 20;
        sets._interpolation_offset    = 0.1;
        sets._interpolation_points    = {400, 10, 200};
        double xi_sth = 1E-3,  eps_sth = 1E-3;
        double xi_pth = 1E-2,  eps_pth = 1E-2;
        double xi_rth = 1E-3,  eps_rth = 1E-3;

        sets._exclusion_offsets   = {1E-1, 3E-1};
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};
        return sets;
    };

    // The P-wave is the dominant isobar
    // In terms of individual isobars this is the only one we need
    class P_wave : public raw_isobar
    {
        // -----------------------------------------------------------------------
        public: 
        
        // Constructor 
        P_wave(isobar_args args) : raw_isobar(args){};

        inline uint    angular_momentum()   { return 1; };
        inline double  phase_shift(double s){ return GKPY::phase_shift(1,1, s); };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            if (iso_id != id::P_wave) return 0.;
            complex  k  = _kinematics->kacser(s), kz = _kinematics->kz(s,t);
            return -3/2*(k*k - kz*kz); // We've multiplied by k^2 which is why singularity_power() = 2
        };
    };
}; // namespace iterateKT 

#endif // PI_ONE_ISOBARS_HPP