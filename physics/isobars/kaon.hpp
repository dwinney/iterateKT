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

#ifndef KAON_ISOBARS_HPP
#define KAON_ISOBARS_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shift.hpp"

namespace iterateKT
{
    // ------------------------------------------------------------------------------
    // All id's for different isobars
 
    // Notation is dIA_IB_C
    // A = two times change in isospin in the decay (1 or 3)
    // B = the total isospin of the 3pi final state (0, 1, 2)
    // C = partial-wave projection of the 2pi subsystem (S0, P1, S2)
    enum class id : unsigned int
    {
                                         //    \delta I    |  total 3pi I  
        dI1_I0_P1,                       //       1/2      |        0      
        dI1_I1_S0, dI1_I1_P1, dI1_I1_S2, //       1/2      |        1 
        dI3_I1_S0, dI3_I1_P1, dI3_I1_S2, //       3/2      |        1 
        dI3_I2_P1, dI3_I2_S2             //       3/2      |        2 
    };

    // ------------------------------------------------------------------------------
    // Isobars

    inline settings default_settings()
    {
        settings sets;
        sets._exclusion_points        = 30;
        sets._exclusion_offsets       = {2.E-2, 2E-2};
        sets._infinitesimal           = 1E-8;
        sets._intermediate_energy     = 1.0;
        sets._cutoff                  = 40.0;
        sets._interpolation_offset    = 1E-4;
        sets._interpolation_points    = {600, 10, 150};

        double xi_sth = 1E-4,   eps_sth = 1E-3;
        double xi_pth = 1E-4,   eps_pth = 1E-3;
        double xi_rth = 1E-3,   eps_rth = 1E-2;
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};
        return sets;
    };

    // Filenames and parameters associated with interpolations of phase_shifts
    phase_args isospin_0 = {"orsay/phase00.dat", 9.9, 1};
    phase_args isospin_1 = {"orsay/phase11.dat", 9.9, 1};
    phase_args isospin_2 = {"orsay/phase02.dat", 9.9, 0};
        
    // ------------------------------------------------------------------------------
    // Isobars

    // \tilde{M}_1
    class I0_P1 : public raw_isobar
    {
        public: 
        I0_P1(isobar_args args) : raw_isobar(args, isospin_1){};

        inline uint    angular_momentum(){ return 1; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI1_I0_P1: return -9*kz*(s-r+kz/3);
                default:            return 0.;
            };
        };
    };

    // M_0 or N_0
    class I1_S0 : public raw_isobar
    {
        public: 
        I1_S0(isobar_args args) : raw_isobar(args, isospin_0){};

        inline uint    angular_momentum(){ return 0; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI1_I1_S0: return 2./3;
                case id::dI1_I1_P1: return 2*(s-r+kz/3);
                case id::dI1_I1_S2: return 20./9;
                default:            return 0;
            };
        };
    };

    // M_1 or N_1
    class I1_P1 : public raw_isobar
    {
        public: 
        I1_P1(isobar_args args) : raw_isobar(args, isospin_1){};

        inline uint    angular_momentum(){ return 1; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI1_I1_S0: return 3*kz;
                case id::dI1_I1_P1: return 9./2*kz*(s-r+kz/3);
                case id::dI1_I1_S2: return -5*kz;
                default:            return 0;
            };
        };
    };

    // M_2 or N_2
    class I1_S2 : public raw_isobar
    {
        public: 
        I1_S2(isobar_args args) : raw_isobar(args, isospin_2){};

        inline uint    angular_momentum(){ return 0; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI1_I1_S0: return 1;
                case id::dI1_I1_P1: return -3./2*(s-r+kz/3);
                case id::dI1_I1_S2: return 1./3;
                default:            return 0;
            };
        };
    };

    // \tilde{N}_1
    class I2_P1 : public raw_isobar
    {
        public: 
        I2_P1(isobar_args args) : raw_isobar(args, isospin_1){};

        inline uint    angular_momentum(){ return 1; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI3_I2_P1: return 9./2*kz*(s-r+kz/3);
                case id::dI3_I2_S2: return -3*kz;
                default:            return 0;
            };
        };
    };

    // \tilde{N}_2
    class I2_S2 : public raw_isobar
    {
        public: 
        I2_S2(isobar_args args) : raw_isobar(args, isospin_2){};

        inline uint    angular_momentum(){ return 0; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::dI3_I2_P1: return -9./2*(s-r+kz/3);
                case id::dI3_I2_S2: return -1;
                default:             return 0;
            };
        };
    };
}; /*  namespace iterateKT */
#endif // KAON_ISOBARS_HPP