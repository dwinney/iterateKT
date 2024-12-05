// Isobars relevant for the decay of isospin-1/2 decay to 3pi in [1]
//
// One may notice that these isobars follow the same structure as those of 
// the eta -> 3pi (e.g. in "isoscalar_pseudoscalar.hpp")
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

#ifndef ISOHALF_PSEUDOSCALAR_HPP
#define ISOHALF_PSEUDOSCALAR_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shift.hpp"

namespace iterateKT
{
    // ------------------------------------------------------------------------------
    // All id's for different isobars
 
    // Notation is dIA_tIB_C
    // A = twice times change in isospin in the decay (1 or 3)
    // B = the total isospin of the 3pi final state (0, 1, 2)
    // C = partial-wave projection of the 2pi subsystem (S0, P1, S2)
    enum class id : unsigned int
    {
                                            //    \delta I    |     total I
        dI1_tI0_P1,                         //       1/2      |        0 
        dI1_tI1_S0, dI1_tI1_P1, dI1_tI1_S2, //       1/2      |        1 
        dI3_tI1_S0, dI3_tI1_P1, dI3_tI1_S2, //       3/2      |        1 
        dI3_tI2_P1, dI3_tI2_S2              //       3/2      |        2 
    };

    // ------------------------------------------------------------------------------
    // Amplitudes

    class charged_mode : public raw_amplitude
    {
        public: 
        charged_mode(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {};
    };

    // ------------------------------------------------------------------------------
    // Isobars

    inline static const settings default_settings()
    {
        settings sets;

        double s0 = M_PION*M_PION;
        sets._exclusion_points        = 10;
        sets._exclusion_offsets       = {2E-2, 3E-2};
        sets._infinitesimal           = 1E-8;
        sets._intermediate_energy     = 1.0;
        sets._cutoff                  = 2.0;
        sets._interpolation_offset    = 1E-4;
        sets._interpolation_points    = {500, 14, 150};

        double xi_sth = 1E-4,     eps_sth = 1E-4;
        double xi_pth = 1E-3,     eps_pth = 1E-3;
        double xi_rth = 0.2*s0,   eps_rth = 1.7*s0;
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};
        return sets;
    };

    // ------------------------------------------------------------------------------
    // Isobars with \delta I = 1/2

    // Total isospin 0
    class dI1_tI0_P1 : public raw_isobar
    {
        public: 
        dI1_tI0_P1(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
        : raw_isobar(xkin, subs, maxsub, sets), _delta1("bern/phase_pipi_1.dat", 80., 1)
        {};

        inline id           get_id()                   { return id::dI1_tI0_P1;  };
        inline unsigned int singularity_power()        { return 2; };
        inline double       phase_shift(double s)      { return _delta1(s/M_PION/M_PION); };
        inline static const settings default_settings(){ return iterateKT::default_settings(); };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI1_tI0_P1: return -3*kz*(3*(s-r)+kz);
                default:             return 0.;
            };
        };
        class phase_shift _delta1;
    };

    // class M0 : public raw_isobar
    // {
    //     public: 
    //     M0(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
    //     : raw_isobar(xkin, subs, maxsub, sets), _delta0("bern/phase_pipi_0.dat", 114., 1)
    //     {};

    //     inline id           get_id()                   { return id::M0;  };
    //     inline unsigned int singularity_power()        { return 0; };
    //     inline double       phase_shift(double s)      { return _delta0(s/M_PION/M_PION); };
    //     inline static const settings default_settings(){ return iterateKT::default_settings(); };
    //     inline complex ksf_kernel(id iso_id, complex s, complex t)
    //     { 
    //         double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
    //         switch (iso_id)
    //         {
    //             case id::M0: return 2/3;
    //             case id::M1: return 2*(s-r+kz/3);
    //             case id::M2: return 20/9;
    //             default:     return 0;
    //         };
    //     };
    //     class phase_shift _delta0;
    // };

    // class M1 : public raw_isobar
    // {
    //     public: 
    //     M1(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
    //     : raw_isobar(xkin, subs, maxsub, sets), _delta1("bern/phase_pipi_1.dat", 80., 1)
    //     {};

    //     inline id           get_id()                   { return id::M1;  };
    //     inline unsigned int singularity_power()        { return 2; };
    //     inline double       phase_shift(double s)      { return _delta1(s/M_PION/M_PION); };
    //     inline static const settings default_settings(){ return iterateKT::default_settings(); };
    //     inline complex ksf_kernel(id iso_id, complex s, complex t)
    //     { 
    //         double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
    //         switch (iso_id)
    //         {
    //             case id::M0: return 3*kz;
    //             case id::M1: return 9/2*kz*(s-r+kz*kz/3);
    //             case id::M2: return -5*kz;
    //             default:     return 0;
    //         };
    //     };
    //     class phase_shift _delta1;
    // };

    // class M2 : public raw_isobar
    // {
    //     public: 
    //     M2(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
    //     : raw_isobar(xkin, subs, maxsub, sets), _delta2("bern/phase_pipi_2.dat", 800., 1)
    //     {};

    //     inline id           get_id()                   { return id::M2;  };
    //     inline unsigned int singularity_power()        { return 0; };
    //     inline double       phase_shift(double s)      { return _delta2(s/M_PION/M_PION); };
    //     inline static const settings default_settings(){ return iterateKT::default_settings(); };
    //     inline complex ksf_kernel(id iso_id, complex s, complex t)
    //     { 
    //         double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
    //         switch (iso_id)
    //         {
    //             case id::M0: return 1;
    //             case id::M1: return -3/2*(s-r+kz/3);
    //             case id::M2: return 1/3;
    //             default:     return 0;
    //         };
    //     };
    //     class phase_shift _delta2;
    // };
}; /*  namespace iterateKT */
#endif // ISOHALF_PSEUDOSCALAR_HPP