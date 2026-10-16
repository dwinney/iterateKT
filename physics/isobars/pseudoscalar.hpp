// Isobars relevant for the decay of isoscalar meson with JP = 0- into 3π
// This allows transitions from I_3π = 0, 1, 2
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

#ifndef PSEUDOSCALAR_ISOBARS_HPP
#define PSEUDOSCALAR_ISOBARS_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shift.hpp"

namespace iterateKT
{
    // ------------------------------------------------------------------------------
    // All id's for different isobars
 
    enum class id : unsigned int // |   I_3π
    {                            // |-----------
        I0_P1,                   // |    0 
        I1_S0, I1_P1, I1_S2,     // |    1
        I2_P1, I2_S2             // |    2
    };

    // ------------------------------------------------------------------------------
    // Isobars with total I_3π = 0

    // I_3π = 0, I_2π = 0, P-wave isobar
    class I0_P1 : public raw_isobar
    {
        public: 
        I0_P1(isobar_args args) : raw_isobar(args){};
        inline unsigned int angular_momentum(){ return 1; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            return (iso_id == id::I0_P1) ? -9*kz*(s-r+kz/3) : 0;
        };
    };

    // ------------------------------------------------------------------------------
    // Isobars with total I_3π = 1
    
    // I_3π = 1, I_2π = 0, S-wave isobar
    class I1_S0 : public raw_isobar
    {
        public: 
        I1_S0(isobar_args args) : raw_isobar(args){};
        inline unsigned int angular_momentum(){ return 0; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::I1_S0: return 2./3;
                case id::I1_P1: return 2*(s-r+kz/3);
                case id::I1_S2: return 20./9;
                default:         return 0;
            };
        };
    };

    // I_3π = 1, I_2π = 1, P-wave
    class I1_P1 : public raw_isobar
    {
        public: 
        I1_P1(isobar_args args) : raw_isobar(args){};
        inline unsigned int angular_momentum(){ return 1; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::I1_S0: return 3*kz;
                case id::I1_P1: return 9./2*kz*(s-r+kz/3);
                case id::I1_S2: return -5*kz;
                default:         return 0;
            };
        };
    };

    // I_3π = 1, I_2π = 2, S-wave
    class I1_S2 : public raw_isobar
    {
        public: 
        I1_S2(isobar_args args) : raw_isobar(args){};
        inline unsigned int angular_momentum(){ return 0; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::I1_S0: return 1;
                case id::I1_P1: return -3./2*(s-r+kz/3);
                case id::I1_S2: return 1./3;
                default:         return 0;
            };
        };
    };

    // ------------------------------------------------------------------------------
    // Isobars with total I_3π = 2

    // I_3π = 2, I_2π = 1, P-wave
    class I2_P1 : public raw_isobar
    {
        public: 
        I2_P1(isobar_args args) : raw_isobar(args){};
        inline unsigned int angular_momentum(){ return 1; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::I2_P1: return 9./2*kz*(s-r+kz/3);
                case id::I2_S2: return -3*kz;
                default:         return 0;
            };
        };
    };

    // I_3π = 2, I_2π = 2, S-wave
    class I2_S2 : public raw_isobar
    {
        public: 
        I2_S2(isobar_args args) : raw_isobar(args){};
        inline unsigned int angular_momentum(){ return 0; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            double  r  = _kinematics->r(); complex kz = _kinematics->kz(s,t);
            switch (iso_id)
            {
                case id::I2_P1: return -9./2*(s-r+kz/3);
                case id::I2_S2: return -1;
                default:         return 0;
            };
        };
    };
}; /*  namespace iterateKT */
#endif // PSEUDOSCALAR_ISOBARS_HPP