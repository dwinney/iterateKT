// Isobars relevant for the decay of meson with JPC = 0-- into 3pi
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

#ifndef PSEUDOSCALAR_HPP
#define PSEUDOSCALAR_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shifts/GKPY.hpp"

namespace iterateKT { namespace pseudoscalar 
{ 
    // Isobar id's
    static const int kI1_S0 = 0;
    static const int kI1_P1 = 1;
    static const int kI1_S2 = 2;

    // Choose default parameters 
    inline static const settings default_settings()
    {
        settings sets;
        sets._adaptive_omnes          = false;
        sets._adaptive_angular        = false;
        sets._adaptive_dispersion     = false;
        sets._infinitesimal           = 1E-5;
        sets._intermediate_energy     = 30;
        sets._cutoff                  = 400;
        sets._interpolation_offset    = 1;
        sets._interpolation_points    = {400, 200, 200};
        double xi_sth = 0.3,  eps_sth = 0.6;
        double xi_pth = 0.9,  eps_pth = 0.9;
        double xi_rth = 0.8,  eps_rth = 1.5;
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};
        return sets;
    };

    // -----------------------------------------------------------------------
    // Isovector \deltaI = 1 amplitude
    // This is isospin violating but C-conserving
    class I1_transition : public raw_amplitude
    {
        public: 
        
        I1_transition(kinematics kin, std::string id) : raw_amplitude(kin,id)
        {};

        inline complex prefactor_s(uint id, complex s, complex t, complex u)
        { 
            switch (id)
            {
                case kI1_S0 : return 1;
                case kI1_S2 : return -2/3;
                default:      return 0;
            };
        };
        inline complex prefactor_t(uint id, complex s, complex t, complex u)
        {
            switch (id)
            {
                case kI1_P1 : return (s-u);
                case kI1_S2 : return 1;
                default:      return 0;
            };
        };
        // Symmetric in t <-> u
        inline complex prefactor_u(uint id, complex s, complex t, complex u){ return prefactor_t(id, s, u, t); };
        
        // Three waves contribute here
        class S0_wave;
        class P1_wave;
        class S2_wave;
    };

    // -----------------------------------------------------------------------
    // I = 0, S-wave isobar
    class I1_transition::S0_wave : public raw_isobar
    {
        public: 
        
        S0_wave(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
        : raw_isobar(xkin, subs, maxsub, sets)
        {};

        inline unsigned int id()                       { return kI1_S0;  };
        inline unsigned int singularity_power()        { return 0; };
        inline double phase_shift(double s)            { return GKPY::phase_shift(0, 0, M_PION*M_PION*s); };
        inline static const settings default_settings(){ return pseudoscalar::default_settings(); };
        inline complex ksf_kernel(int iso_id, complex s, complex t)
        { 
            double  m2 = _kinematics->m2(), M2 = _kinematics->M2(), r  = M2/3 + m2;
            complex k  = _kinematics->kacser(s), kz = 2*t+s-M2-3*m2;
            switch (iso_id)
            {
                case kI1_S0: return 2./3;
                case kI1_P1: return 2.*((s-r)+kz/3);
                case kI1_S2: return 20./9;
                default:     return 0.;
            };
        };
    };

    // -----------------------------------------------------------------------
    // I = 1, P-wave isobar
    class I1_transition::P1_wave : public raw_isobar
    {
        public: 
        
        P1_wave(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
        : raw_isobar(xkin, subs, maxsub, sets)
        {};

        inline unsigned int id()                       { return kI1_P1;  };
        inline unsigned int singularity_power()        { return 2; };
        inline double       phase_shift(double s)      { return GKPY::phase_shift(1, 1, M_PION*M_PION*s); };
        inline static const settings default_settings(){ return pseudoscalar::default_settings(); };

        inline complex ksf_kernel(int iso_id, complex s, complex t)
        { 
            double  m2 = _kinematics->m2(), M2 = _kinematics->M2(), r  = M2/3 + m2;
            complex k  = _kinematics->kacser(s), kz = 2*t+s-M2-3*m2;
            switch (iso_id)
            {
                case kI1_S0: return 3*kz;
                case kI1_P1: return 3*kz/2*(3*(s-r)+kz);
                case kI1_S2: return -5*kz;
                default:     return 0.;
            };
        };
    };

    // -----------------------------------------------------------------------
    // I = 2, S-wave isobar
    class I1_transition::S2_wave : public raw_isobar
    {
        public: 
        
        S2_wave(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
        : raw_isobar(xkin, subs, maxsub, sets)
        {};

        inline unsigned int id()                       { return kI1_S2;  };
        inline unsigned int singularity_power()        { return 0; };
        inline double phase_shift(double s)            { return GKPY::phase_shift(2, 0, M_PION*M_PION*s); };
        inline static const settings default_settings(){ return pseudoscalar::default_settings(); };
        inline complex ksf_kernel(int iso_id, complex s, complex t)
        { 
            double  m2 = _kinematics->m2(), M2 = _kinematics->M2(), r  = M2/3 + m2;
            complex k  = _kinematics->kacser(s), kz = 2*t+s-M2-3*m2;
            switch (iso_id)
            {
                case kI1_S0: return 1.;
                case kI1_P1: return -(3*(s-r)+kz)/2;
                case kI1_S2: return 2./3;
                default:    return 0.;
            };
        };
    };

}; /*  namespace iterateKT */ }; // namespace pseudoscalar

#endif // PSEUDOSCALAR_HPP