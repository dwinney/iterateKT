// Amplitudes relevant for the decay of isoscalar meson with JP = 0- into 3pi
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] -  https://arxiv.org/abs/2111.02417
// ------------------------------------------------------------------------------

#ifndef ETA_AMPLITUDES_HPP
#define ETA_AMPLITUDES_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shift.hpp"

namespace iterateKT
{
    // Isospin breaking normalization
    constexpr double XI = -0.140;

    class charged_mode : public raw_amplitude
    {
        public: 
        charged_mode(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {};
        
        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            switch(iso_id)
            {
                case (id::dI0_P1): return (t-u);
                case (id::dI1_S0): return XI;
                case (id::dI1_S2): return XI*(-2/3);
                case (id::dI2_P1): return 2*(u-t);
                default: return 0;
            };
        };

        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            switch(iso_id)
            {
                case (id::dI0_P1): return (u-s);
                case (id::dI1_P1): return XI*(s-u);
                case (id::dI1_S2): return XI;
                case (id::dI2_P1): return (u-s);
                case (id::dI2_S2): return -1;
                default: return 0;
            };
        };

        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            switch(iso_id)
            {
                case (id::dI0_P1): return (s-t);
                case (id::dI1_P1): return XI*(s-t);
                case (id::dI1_S2): return XI;
                case (id::dI2_P1): return (s-t);
                case (id::dI2_S2): return 1;
                default: return 0;
            };
        };
    };

    class neutral_mode : public raw_amplitude
    {
        public: 
        neutral_mode(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {};
        
        // 3 identical particles in the final state
        inline combinatorial_factor(){ return 6; };

        // Only S-waves contribute by Bose symmetry.
        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            switch(iso_id)
            {
                case (id::dI1_S0): return XI;
                case (id::dI1_S2): return XI*4/3;
                default: return 0;
            };
        };
        // The rest are symmetric
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u){ return prefactor_t(iso_id, t, s, u); };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u){ return prefactor_t(iso_id, u, t, s); };
    };
};

#endif