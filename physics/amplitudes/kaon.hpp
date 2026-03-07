// Isobars relevant for the decay of isospin-1/2 pseudoscalar decay to 3π
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2403.17570
// [2] - https://arxiv.org/abs/2111.02417
// ------------------------------------------------------------------------------

#ifndef KAON_AMPLITUDES_HPP
#define KAON_AMPLITUDES_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shift.hpp"
#include "amplitude.hpp"
#include "isobars/pseudoscalar.hpp"

// For a general K_Pi1Pi2Pi3, we have:
// s = (K - π1)^2 = (π2 + π3)^2 = s1
// t = (K - π2)^2 = (π1 + π3)^2 = s2
// u = (K - π3)^2 = (π1 + π2)^2 = s3

// We assume the symmetric channels are t <-> u

namespace iterateKT
{
    inline settings default_settings()
    {
        settings sets;
        sets._derivative_h            = 1E-6;
        sets._exclusion_points        = 20;
        sets._exclusion_offsets       = {2.E-2, 5E-2};
        sets._infinitesimal           = 1E-8;
        sets._intermediate_energy     = 1.0;
        sets._cutoff                  = 20.0;
        sets._interpolation_offset    = 1E-4;
        sets._interpolation_points    = {400, 10, 300};
        sets._extra_cusp              = norm(0.980);

        double xi_sth = 1E-3,   eps_sth = 1E-3;
        double xi_pth = 2E-3,   eps_pth = 1E-2;
        double xi_rth = 2E-2,   eps_rth = 2E-2;
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};

        std::string path = main_dir() + "/physics/phase_shifts/";
        phase_args iso_0 = {path+"madrid/delta_00.dat", 1.69,  1, 1, 2};
        phase_args iso_1 = {path+"madrid/delta_11.dat", 1.69,  1, 1, 2};
        phase_args iso_2 = {path+"madrid/delta_02.dat", 9.99,  0, 1, 2};
                
        sets._phase_shifts = { {id::I0_P1, iso_1}, 
                               {id::I1_S0, iso_0}, {id::I1_P1, iso_1}, {id::I1_S2, iso_2},
                               {id::I2_P1, iso_1}, {id::I2_S2, iso_2}};
        return sets;
    };
    
    // ------------------------------------------------------------------------------
    // These are the invariant amplitudes from a 3π state of total isospin

    // I_3π = 1 amplitude
    class I1 : public raw_amplitude
    {
        public:

        I1(kinematics kin) : raw_amplitude(kin) {};

        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case id::I1_S0: return 1;
                case id::I1_S2: return -2./3;
                default: return 0;
            };
        };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case id::I1_P1: return (s-u);
                case id::I1_S2: return 1;
                default: return 0;
            };
        };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            return prefactor_t(iso_id, s, u, t);
        };
    };
    
    // I_3π = 2 amplitude
    class I2 : public raw_amplitude
    {
        public:

        I2(kinematics kin) : raw_amplitude(kin) {};

        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            return (iso_id == id::I2_S2) ? 1 : 0;
        };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case id::I2_P1: return -3*(s-u)/2;
                case id::I2_S2: return -1./2;
                default: return 0;
            };
        };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            return prefactor_t(iso_id, s, u, t);
        };
    };

    // ------------------------------------------------------------------------------
    // These are the physical amplitudes in the charged basis
    
    // Allowed to access both K+ -> pi+pi+pi- and pi0pi0pi+
    enum class option : unsigned int { P_ppm, P_zzp };

    class charged_kaon : public raw_amplitude
    {
        public: 
        
        charged_kaon(kinematics xkin) : raw_amplitude(xkin)
        {
            _F = new_amplitude<I1>(xkin); _H = new_amplitude<I2>(xkin);
        };

        inline void set_option(option opt){ _charged = (opt == option::P_ppm); };

        // P_ppm(s,t,u) = F(t,s,u) + F(u,t,s) + H(s,t,u)
        // P_zzp(s,t,u) = F(s,t,u)            + H(s,t,u)

        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            complex F = (_charged) ?  _F->prefactor_t(iso_id, t, s, u) + _F->prefactor_u(iso_id, u, t, s)
                                   :  _F->prefactor_s(iso_id, s, t, u);
            complex H = _H->prefactor_s(iso_id, s, t, u);

            return F + H;
        };

        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            complex F = (_charged) ?  _F->prefactor_s(iso_id, t, s, u) + _F->prefactor_t(iso_id, u, t, s)
                                   :  _F->prefactor_t(iso_id, s, t, u);
            complex H = _H->prefactor_t(iso_id, s, t, u);

            return F + H;
        };

        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            complex F = (_charged) ?  _F->prefactor_s(iso_id, u, t, s) + _F->prefactor_u(iso_id, t, s, u)
                                   :  _F->prefactor_u(iso_id, s, t, u);
            complex H = _H->prefactor_u(iso_id, s, t, u);

            return F + H;
        };

        // Two identical particles
        inline double combinatorial_factor(){ return 2; }; 

        private: 
        amplitude _F, _H;
        bool _charged = true;
    };
};

#endif