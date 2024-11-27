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

#ifndef PSEUDOSCALAR_HPP
#define PSEUDOSCALAR_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shift.hpp"

namespace iterateKT { namespace pseudoscalar
{
    // Isobar id's
    static const int kdI0_P1 = 0;
    static const int kdI1_S0 = 1; 
    static const int kdI1_P1 = 2;
    static const int kdI1_S2 = 3;
    static const int kdI2_P1 = 4;
    static const int kdI2_S2 = 5;

    inline static const settings default_settings()
    {
        settings sets;
        sets._adaptive_omnes          = false;
        sets._adaptive_angular        = false;
        sets._adaptive_cauchy         = true; 
        sets._adaptive_pseudo         = true; 

        sets._pseudo_depth = 10;
        sets._cauchy_depth = 5;

        sets._interpolation_type = ROOT::Math::Interpolation::Type::kCSPLINE;

        sets._exclusion_points        = 10;
        sets._exclusion_offsets       = {1, 3};
        sets._infinitesimal           = 1E-5;
        sets._intermediate_energy     = 100;
        sets._cutoff                  = 1000;
        sets._interpolation_offset    = 0.1;
        sets._interpolation_points    = {300, 16, 100};

        double xi_sth = 0.3,   eps_sth = 0.3;
        double xi_pth = 0.5,   eps_pth = 0.5;
        double xi_rth = 1.5,   eps_rth = 1.5;
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};
        return sets;
    };

    // ------------------------------------------------------------------------------
    // Amplitudes

    class charged_mode : public raw_amplitude
    {
        public: 
        charged_mode(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {};
    };

    class neutral_mode : public raw_amplitude
    {
        public: 
        neutral_mode(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {};
    };


    // ------------------------------------------------------------------------------
    // Isobars

    // dI = 0, I = 0, P-wave isobar
    class dI0_P1 : public raw_isobar
    {
        public: 
        dI0_P1(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
        : raw_isobar(xkin, subs, maxsub, sets), _delta1("bern/phase_pipi_1.dat", 80., 1)
        {};

        inline unsigned int id()                       { return kdI0_P1;  };
        inline unsigned int singularity_power()        { return 2; };
        inline double       phase_shift(double s)      { return _delta1(s); };
        inline static const settings default_settings(){ return pseudoscalar::default_settings(); };
        inline complex ksf_kernel(int iso_id, complex s, complex t)
        { 
            double  m2 = _kinematics->m2(), M2 = _kinematics->M2(), r  = M2/3 + m2;
            complex k  = _kinematics->kacser(s), kz = 2*t+s-M2-3*m2;
            switch (iso_id)
            {
                case kdI0_P1: return -3*kz*(3*(s-r)+kz);
                default:      return 0.;
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

        inline unsigned int id()                       { return kdI1_S0;  };
        inline unsigned int singularity_power()        { return 0; };
        inline double phase_shift(double s)            { return _delta0(s); };
        inline static const settings default_settings(){ return pseudoscalar::default_settings(); };
        inline complex ksf_kernel(int iso_id, complex s, complex t)
        { 
            double  m2 = _kinematics->m2(), M2 = _kinematics->M2(), r  = M2/3 + m2;
            complex k  = _kinematics->kacser(s), kz = 2*t+s-M2-3*m2;
            switch (iso_id)
            {
                case kdI1_S0: return 2./3;
                case kdI1_P1: return 2.*((s-r)+kz/3);
                case kdI1_S2: return 20./9;
                default:      return 0.;
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

        inline unsigned int id()                       { return kdI1_P1;  };
        inline unsigned int singularity_power()        { return 2; };
        inline double       phase_shift(double s)      { return _delta1(s); };
        inline static const settings default_settings(){ return pseudoscalar::default_settings(); };
        inline complex ksf_kernel(int iso_id, complex s, complex t)
        { 
            double  m2 = _kinematics->m2(), M2 = _kinematics->M2(), r  = M2/3 + m2;
            complex k  = _kinematics->kacser(s), kz = 2*t+s-M2-3*m2;
            switch (iso_id)
            {
                case kdI1_S0: return 3*kz;
                case kdI1_P1: return 3*kz/2*(3*(s-r)+kz);
                case kdI1_S2: return -5*kz;
                default:      return 0.;
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

        inline unsigned int id()                       { return kdI1_S2;  };
        inline unsigned int singularity_power()        { return 0; };
        inline double phase_shift(double s)            { return _delta2(s); };
        inline static const settings default_settings(){ return pseudoscalar::default_settings(); };
        inline complex ksf_kernel(int iso_id, complex s, complex t)
        { 
            double  m2 = _kinematics->m2(), M2 = _kinematics->M2(), r  = M2/3 + m2;
            complex k  = _kinematics->kacser(s), kz = 2*t+s-M2-3*m2;
            switch (iso_id)
            {
                case kdI1_S0: return 1.;
                case kdI1_P1: return -(3*(s-r)+kz)/2;
                case kdI1_S2: return 2./3;
                default:      return 0.;
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

        inline unsigned int id()                       { return kdI2_P1;  };
        inline unsigned int singularity_power()        { return 2; };
        inline double       phase_shift(double s)      { return _delta1(s); };
        inline static const settings default_settings(){ return pseudoscalar::default_settings(); };
        inline complex ksf_kernel(int iso_id, complex s, complex t)
        { 
            double  m2 = _kinematics->m2(), M2 = _kinematics->M2(), r  = M2/3 + m2;
            complex k  = _kinematics->kacser(s), kz = 2*t+s-M2-3*m2;
            switch (iso_id)
            {
                case kdI2_P1: return 3*kz/2*(3*(s-r)+kz);
                case kdI2_S2: return 3*kz;
                default:      return 0.;
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

        inline unsigned int id()                       { return kdI2_S2;  };
        inline unsigned int singularity_power()        { return 0; };
        inline double       phase_shift(double s)      { return _delta2(s); };
        inline static const settings default_settings(){ return pseudoscalar::default_settings(); };
        inline complex ksf_kernel(int iso_id, complex s, complex t)
        { 
            double  m2 = _kinematics->m2(), M2 = _kinematics->M2(), r  = M2/3 + m2;
            complex k  = _kinematics->kacser(s), kz = 2*t+s-M2-3*m2;
            switch (iso_id)
            {
                case kdI2_P1: return 3./2*(3*(s-r)+kz);
                case kdI2_S2: return -1.;
                default:      return 0.;
            };
        };
        class phase_shift _delta2;
    };

}; /*  namespace iterateKT */ }; // namespace pseudoscalar
#endif // PSEUDOSCALAR_HPP