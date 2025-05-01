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

#ifndef ETA_ISOBARS_HPP
#define ETA_ISOBARS_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shift.hpp"

namespace iterateKT
{
    // ------------------------------------------------------------------------------
    // All id's for different isobars
 
    // All isobars couple the 3pi state to the isospin-0
    // Since all particles are iso-bosons, we only need to specify the \delta I
    enum class id : unsigned int
    {   
                                   //  Isospin violating?  |  C-violating?
        dI0_P1,                    //        no            |      yes
        dI1_S0, dI1_P1, dI1_S2,    //        yes           |      no
        dI2_P1, dI2_S2             //        yes           |      yes
    };

    // All our default settings related to technicals of KT solution
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
        return sets;
    };

    // Filenames and parameters associated with interpolations of phase_shifts
    phase_args isospin_0 = {"bern/phase_pipi_0.dat", 114., 1};
    phase_args isospin_1 = {"bern/phase_pipi_1.dat", 80,   1};
    phase_args isospin_2 = {"bern/phase_pipi_2.dat", 800,  0};

    // ------------------------------------------------------------------------------
    // Isobars with total \delta I = 0
    // These are isospin conserving but C-violating

    // dI = 0, I = 0, P-wave isobar
    class dI0_P1 : public raw_isobar
    {
        public: 
        dI0_P1(isobar_args args) : raw_isobar(args, isospin_1){};
        inline unsigned int angular_momentum(){ return 1; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI0_P1: return -9*kz*(s-r+kz/3);
                default:         return 0;
            };
        };
    };

    // ------------------------------------------------------------------------------
    // Isobars with total \delta I = 1
    // These are C-conserving but isospin violating
    
    // dI = 1, I = 0, S-wave isobar
    class dI1_S0 : public raw_isobar
    {
        public: 
        dI1_S0(isobar_args args) : raw_isobar(args, isospin_0){};
        inline unsigned int angular_momentum(){ return 0; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI1_S0: return 2./3;
                case id::dI1_P1: return 2*(s-r+kz/3);
                case id::dI1_S2: return 20./9;
                default:         return 0;
            };
        };
        class phase_shift _delta0;
    };

    // dI = 1, I = 1, P-wave
    class dI1_P1 : public raw_isobar
    {
        public: 
        dI1_P1(isobar_args args) : raw_isobar(args, isospin_1){};
        inline unsigned int angular_momentum(){ return 1; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI1_S0: return 3*kz;
                case id::dI1_P1: return 9./2*kz*(s-r+kz/3);
                case id::dI1_S2: return -5*kz;
                default:         return 0;
            };
        };
    };

    // dI = 1, I = 2, S-wave
    class dI1_S2 : public raw_isobar
    {
        public: 
        dI1_S2(isobar_args args) : raw_isobar(args, isospin_2){};
        inline unsigned int angular_momentum(){ return 0; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI1_S0: return 1;
                case id::dI1_P1: return -3./2*(s-r+kz/3);
                case id::dI1_S2: return 1./3;
                default:         return 0;
            };
        };
    };

    // ------------------------------------------------------------------------------
    // Isobars with total \delta I = 2
    // These are both C and isospin violating

    // dI = 2, I = 1, P-wave
    class dI2_P1 : public raw_isobar
    {
        public: 
        dI2_P1(isobar_args args) : raw_isobar(args, isospin_1){};
        inline unsigned int angular_momentum(){ return 1; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI2_P1: return 9./2*kz*(s-r+kz/3);
                case id::dI2_S2: return 3*kz;
                default:         return 0;
            };
        };
        class phase_shift _delta1;
    };

    // dI = 2, I = 2, S-wave
    class dI2_S2 : public raw_isobar
    {
        public: 
        dI2_S2(isobar_args args) : raw_isobar(args, isospin_2){};
        inline unsigned int angular_momentum(){ return 0; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI2_P1: return 9./2*(s-r+kz/3);
                case id::dI2_S2: return -1;
                default:         return 0;
            };
        };
    };
}; /*  namespace iterateKT */
#endif // ETA_ISOBARS_HPP