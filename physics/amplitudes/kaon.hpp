// Isobars relevant for the decay of isospin-1/2 decay to 3pi in [1]
//
// One may notice that these isobars follow the same structure as those of 
// the eta -> 3pi (e.g. in "eta.hpp")
// but these are redefined here anyway to avoid confusion
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2403.17570
// ------------------------------------------------------------------------------

#ifndef KAON_AMPLITUDES_HPP
#define KAON_AMPLITUDES_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shift.hpp"
#include "isobars/kaon.hpp"

// The amplitudes are named with respect to isospin projections of the decay particle
// into three pions. The order matters in that the definitions of s, t, and u.

// For a general K_Pi1Pi2Pi3, we have:
// s = (K - Pi1)^2 = (Pi2 + Pi3)^2
// t = (K - Pi2)^2 = (Pi1 + Pi3)^2
// u = (K - Pi3)^2 = (Pi1 + Pi2)^2

namespace iterateKT
{
    //--------------------------------------------------------------------------
    // K+ -> pi+ pi+ pi-         
    class Kp_PipPipPim : public raw_amplitude
    {
        public: 
        Kp_PipPipPim(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {};

        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_S0): return -1;
                case (id::dI1_I1_P1): return -(t-u);
                case (id::dI1_I1_S2): return -1/3;
                case (id::dI3_I1_S0): return -1;
                case (id::dI3_I1_P1): return -(t-u);
                case (id::dI3_I1_S2): return -1/3;
                case (id::dI3_I2_P1): return -3/2*(t-u);
                case (id::dI3_I2_S2): return +1/2;
                default: return 0;
            };
        };

        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            return prefactor_s(iso_id, t, s, u);
        };

        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_S2): return -2;
                case (id::dI3_I1_S2): return -2;
                case (id::dI3_I2_S2): return -1;
                default: return 0;
            };
        };
    };

    //--------------------------------------------------------------------------
    // K+ -> pi0 pi0 pi+         
    class Kp_PizPizPip : public raw_amplitude
    {
        public: 
        Kp_PizPizPip(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {};

        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_P1): return -    (t-u);
                case (id::dI1_I1_S2): return +1;
                case (id::dI3_I1_P1): return -    (t-u);
                case (id::dI3_I1_S2): return +1;
                case (id::dI3_I2_P1): return +3/2*(t-u);
                case (id::dI3_I2_S2): return -1/2;
                default: return 0;
            };
        };

        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            return prefactor_s(iso_id, t, s, u);
        };

        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_S0): return +1;
                case (id::dI1_I1_S2): return -2/3;
                case (id::dI3_I1_S0): return +1;
                case (id::dI3_I1_S2): return -2/3;
                case (id::dI3_I2_S2): return +1;
                default: return 0;
            };
        };
    };

    
    //--------------------------------------------------------------------------
    // KL -> pi+ pi- pi0         
    class KL_PipPimPiz : public raw_amplitude
    {
        public: 
        KL_PipPimPiz(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {};

        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_P1): return +  (t-u);
                case (id::dI1_I1_S2): return -1;
                case (id::dI3_I1_P1): return -2*(t-u);
                case (id::dI3_I1_S2): return +2;
                default: return 0;
            };
        };

        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            return prefactor_s(iso_id, t, s, u);
        };

        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_S0): return -1;
                case (id::dI1_I1_S2): return +2/3;
                case (id::dI3_I1_S0): return +2;
                case (id::dI3_I1_S2): return -4/3;
                default: return 0;
            };
        };
    };

    //--------------------------------------------------------------------------
    // KS -> pi+ pi- pi0         
    class KS_PipPimPiz : public raw_amplitude
    {
        public: 
        KS_PipPimPiz(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {};

        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case (id::dI1_I0_P1): return +(t-u);
                case (id::dI3_I2_P1): return +(t-u);
                case (id::dI3_I2_S2): return +1;
                default: return 0;
            };
        };

        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            return -prefactor_s(iso_id, t, s, u);
        };

        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case (id::dI1_I0_P1): return +  (s-t);
                case (id::dI3_I2_P1): return -2*(s-t);
                default: return 0;
            };
        };
    };

    
    //--------------------------------------------------------------------------
    // KL -> pi0 pi0 pi0         
    class KL_PizPizPiz : public raw_amplitude
    {
        public: 
        KL_PizPizPiz(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {};

        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_S0): return +1;
                case (id::dI1_I1_S2): return +4/3;
                case (id::dI3_I1_S0): return -2;
                case (id::dI3_I1_S2): return -8/3;
                default: return 0;
            };
        };

        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            return prefactor_s(iso_id, t, s, u);
        };

        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            return prefactor_s(iso_id, u, t, s);
        };
    };

    //--------------------------------------------------------------------------
    // Generic K -> 3 pi this will contain all the above which can be accessed via the set_option() function
    
    enum class option : unsigned int
    {
        Kp_PipPipPim, Kp_PizPizPip, KL_PipPimPiz, KS_PipPimPiz, KL_PizPizPiz
    };

    class K_3Pi : public raw_amplitude
    {
        public:
        K_3Pi(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {
            p_ppm = new_amplitude<Kp_PipPipPim>(xkin);
            p_zzp = new_amplitude<Kp_PizPizPip>(xkin);
            L_pmz = new_amplitude<KL_PipPimPiz>(xkin);
            S_pmz = new_amplitude<KS_PipPimPiz>(xkin);
            L_zzz = new_amplitude<KL_PizPizPiz>(xkin);
            current = p_ppm;
        };

        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            return current->prefactor_s(iso_id, s, t, u);
        };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            return current->prefactor_t(iso_id, s, t, u);
        };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            return current->prefactor_u(iso_id, s, t, u);
        };

        inline void set_option(option opt)
        {
            switch (opt)
            {
                case (option::Kp_PipPipPim): current = p_ppm; return;
                case (option::Kp_PizPizPip): current = p_zzp; return;
                case (option::KL_PipPimPiz): current = L_pmz; return;
                case (option::KS_PipPimPiz): current = S_pmz; return;
                case (option::KL_PizPizPiz): current = L_zzz; return;
                default: return;
            };
        };

        private:

        amplitude current; 
        amplitude p_ppm, p_zzp, L_pmz, S_pmz, L_zzz;
    };
};

#endif