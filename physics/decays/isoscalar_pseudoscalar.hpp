// Isobars relevant for the decay of isoscalar meson with JP = 0- into 3pi
// This allows delta I = 0, 1, 2 transitions
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

#ifndef ISOSCALAR_PSEUDOSCALAR_HPP
#define ISOSCALAR_PSEUDOSCALAR_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shift.hpp"

namespace iterateKT
{
    // ------------------------------------------------------------------------------
    // All id's for different isobars
 
    enum class id : unsigned int
    {
        dI0_P1,
        dI1_S0, dI1_P1, dI1_S2,
        dI2_P1, dI2_S2
    };

    // ------------------------------------------------------------------------------
    // Amplitudes


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

    // ------------------------------------------------------------------------------
    // Isobars involved in the transitions with integer total \delta I 

    inline static const settings default_settings()
    {
        settings sets;
        
        // All depths to zero, we use non-adaptive integration
        sets._pseudo_integrator_depth  = 0;
        sets._cauchy_integrator_depth  = 0;
        sets._angular_integrator_depth = 0;
        sets._omnes_integrator_depth   = 0;

        sets._interpolation_type = ROOT::Math::Interpolation::Type::kCSPLINE;

        sets._exclusion_points        = 10;
        sets._exclusion_offsets       = {2, 3};
        sets._infinitesimal           = 1E-5;
        sets._intermediate_energy     = 100;
        sets._cutoff                  = 1000;
        sets._interpolation_offset    = 0.1;
        sets._interpolation_points    = {300, 22, 150};

        double xi_sth = 0.2,   eps_sth = 0.2;
        double xi_pth = 0.5,   eps_pth = 0.5;
        double xi_rth = 0.2,   eps_rth = 1.7;
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};
        return sets;
    };

    // dI = 0, I = 0, P-wave isobar
    class dI0_P1 : public raw_isobar
    {
        public: 
        dI0_P1(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
        : raw_isobar(xkin, subs, maxsub, sets), _delta1("bern/phase_pipi_1.dat", 80., 1)
        {};

        inline id           get_id()                   { return id::dI0_P1;  };
        inline unsigned int singularity_power()        { return 2; };
        inline double       phase_shift(double s)      { return _delta1(s); };
        inline static const settings default_settings(){ return iterateKT::default_settings(); };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI0_P1: return -3*kz*(3*(s-r)+kz);
                default:         return 0.;
            };
        };
        class phase_shift _delta1;
    };

    // dI = 1, I = 0, S-wave isobar
    class dI1_S0 : public raw_isobar
    {
        public: 
        dI1_S0(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
        : raw_isobar(xkin, subs, maxsub, sets), _delta0("bern/phase_pipi_0.dat", 114., 1)
        {};

        inline id           get_id()                   { return id::dI1_S0;  };
        inline unsigned int singularity_power()        { return 0; };
        inline double phase_shift(double s)            { return _delta0(s); };
        inline static const settings default_settings(){ return iterateKT::default_settings(); };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI1_S0: return 2./3;
                case id::dI1_P1: return 2.*((s-r)+kz/3);
                case id::dI1_S2: return 20./9;
                default:         return 0.;
            };
        };
        class phase_shift _delta0;
    };

    // dI = 1, I = 1, P-wave
    class dI1_P1 : public raw_isobar
    {
        public: 
        dI1_P1(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
        : raw_isobar(xkin, subs, maxsub, sets), _delta1("bern/phase_pipi_1.dat", 80., 1)
        {};

        inline id           get_id()                   { return id::dI1_P1;  };
        inline unsigned int singularity_power()        { return 2; };
        inline double       phase_shift(double s)      { return _delta1(s); };
        inline static const settings default_settings(){ return iterateKT::default_settings(); };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI1_S0: return 3*kz;
                case id::dI1_P1: return 3*kz/2*(3*(s-r)+kz);
                case id::dI1_S2: return -5*kz;
                default:         return 0.;
            };
        };
        class phase_shift _delta1;
    };

    // dI = 1, I = 2, S-wave
    class dI1_S2 : public raw_isobar
    {
        public: 
        dI1_S2(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
        : raw_isobar(xkin, subs, maxsub, sets), _delta2("bern/phase_pipi_2.dat", 800, 0)
        {};

        inline id           get_id()                   { return id::dI1_S2;  };
        inline unsigned int singularity_power()        { return 0; };
        inline double phase_shift(double s)            { return _delta2(s); };
        inline static const settings default_settings(){ return iterateKT::default_settings(); };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI1_S0: return 1.;
                case id::dI1_P1: return -(3*(s-r)+kz)/2;
                case id::dI1_S2: return 2./3;
                default:         return 0.;
            };
        };
        class phase_shift _delta2;
    };

    // dI = 2, I = 1, P-wave
    class dI2_P1 : public raw_isobar
    {
        public: 
        dI2_P1(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
        : raw_isobar(xkin, subs, maxsub, sets), _delta1("bern/phase_pipi_1.dat", 80., 1)
        {};

        inline id           get_id()                   { return id::dI2_P1;  };
        inline unsigned int singularity_power()        { return 2; };
        inline double       phase_shift(double s)      { return _delta1(s); };
        inline static const settings default_settings(){ return iterateKT::default_settings(); };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI2_P1: return 3*kz/2*(3*(s-r)+kz);
                case id::dI2_S2: return 3*kz;
                default:         return 0.;
            };
        };
        class phase_shift _delta1;
    };

    // dI = 2, I = 2, S-wave
    class dI2_S2 : public raw_isobar
    {
        public: 
        dI2_S2(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
        : raw_isobar(xkin, subs, maxsub, sets), _delta2("bern/phase_pipi_2.dat", 800, 0)
        {};

        inline id           get_id()                   { return id::dI2_S2;  };
        inline unsigned int singularity_power()        { return 0; };
        inline double       phase_shift(double s)      { return _delta2(s); };
        inline static const settings default_settings(){ return iterateKT::default_settings(); };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI2_P1: return 3./2*(3*(s-r)+kz);
                case id::dI2_S2: return -1.;
                default:         return 0.;
            };
        };
        class phase_shift _delta2;
    };

}; /*  namespace iterateKT */
#endif // ETA_HPP