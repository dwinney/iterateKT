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
#include "isobars/pseudoscalar.hpp"

namespace iterateKT
{
    // Isospin breaking normalization
    constexpr double XI = -0.140;

    // Default settings for isobars
    inline settings default_settings()
    {
        settings sets;
        sets._exclusion_points        = 10;
        sets._exclusion_offsets       = {2, 3};
        sets._infinitesimal           = 1E-5;
        sets._intermediate_energy     = 100;
        sets._cutoff                  = 1000;
        sets._interpolation_offset    = 0.1;
        sets._interpolation_points    = {400, 22, 150};
        double xi_sth = 0.4,   eps_sth = 0.4;
        double xi_pth = 0.4,   eps_pth = 0.4;
        double xi_rth = 0.3,   eps_rth = 1.9;
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};
        
        phase_args iso_0 = {"bern/phase_pipi_0.dat", 114., 1, 1, 2};
        phase_args iso_1 = {"bern/phase_pipi_1.dat", 800., 1, 1, 2};
        phase_args iso_2 = {"bern/phase_pipi_2.dat", 800,  0, 1, 2};
        sets._phase_shifts = { {id::I0_P1, iso_1}, 
                               {id::I1_S0, iso_0}, {id::I1_P1, iso_1}, {id::I1_S2, iso_2},
                               {id::I2_P1, iso_1}, {id::I2_S2, iso_2}};
        return sets;
    };

    class charged_mode : public raw_amplitude
    {
        public: 
        charged_mode(kinematics xkin) : raw_amplitude(xkin) {};
        
        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            switch(iso_id)
            {
                case (id::I0_P1): return (t-u);
                case (id::I1_S0): return XI;
                case (id::I1_S2): return XI*(-2/3);
                case (id::I2_P1): return 2*(u-t);
                default: return 0;
            };
        };

        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            switch(iso_id)
            {
                case (id::I0_P1): return (u-s);
                case (id::I1_P1): return XI*(s-u);
                case (id::I1_S2): return XI;
                case (id::I2_P1): return (u-s);
                case (id::I2_S2): return -1;
                default: return 0;
            };
        };

        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            switch(iso_id)
            {
                case (id::I0_P1): return (s-t);
                case (id::I1_P1): return XI*(s-t);
                case (id::I1_S2): return XI;
                case (id::I2_P1): return (s-t);
                case (id::I2_S2): return 1;
                default: return 0;
            };
        };
    };

    class neutral_mode : public raw_amplitude
    {
        public: 
        neutral_mode(kinematics xkin) : raw_amplitude(xkin) {};
        
        // 3 identical particles in the final state
        inline double combinatorial_factor(){ return 6; };

        // Only S-waves contribute by Bose symmetry.
        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            switch(iso_id)
            {
                case (id::I1_S0): return XI;
                case (id::I1_S2): return XI*4/3;
                default: return 0;
            };
        };
        // The rest are symmetric
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u){ return prefactor_t(iso_id, t, s, u); };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u){ return prefactor_t(iso_id, u, t, s); };
    };
};

#endif