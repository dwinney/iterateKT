// Isobars relevant for the decay of meson with JPC = 1-- into 3pi as in Ref. [1]
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

#ifndef V_TO_3PI_HPP
#define V_TO_3PI_HPP

#include "isobar.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phaseshifts/GKPY.hpp"

namespace iterateKT { namespace V_to_3pi
{
    // Isobar id's
    static const int kP_wave = 0;

    // This defines the full amplitude, i.e. how the isobars are combined
    // Here is where we usually put the isospin combinations etc
    class isoscalar : public raw_amplitude
    {
        public: 
        
        // Constructor
        isoscalar(kinematics kin, std::string id) : raw_amplitude(kin,id)
        {};

        // We have no kinematic factors and only one isobar so simply return 1.
        inline complex s_channel_prefactor(uint id, complex s, complex t, complex u)
        {
            return (id == kP_wave) ? 1. : 0.;
        };

        // We're completely symmetric here so these are the same
        inline complex t_channel_prefactor(uint id, complex s, complex t, complex u)
        {
            return (id == kP_wave) ? 1. : 0.;
        };
        inline complex u_channel_prefactor(uint id, complex s, complex t, complex u)
        {
            return (id == kP_wave) ? 1. : 0.;
        };
    };
    
    // The P-wave is the dominant isobar
    // In terms of individual isobars this is the only one we need
    class P_wave : public raw_isobar
    {
        // -----------------------------------------------------------------------
        public: 
        
        // Constructor 
        P_wave(kinematics xkin, int nsub, settings sets) : raw_isobar(xkin, nsub, sets)
        {};

        // Virtual functions
        inline unsigned int  id() { return kP_wave;  };
        inline std::string name() { return "P-wave"; };

        // Because the P-wave involes a sintheta = 1-z^2, we have two power of 1/kappa
        // which lead to pseudo threshold singularities
        // The TOTAL singularity power is always +1 from this (one factor from the jacobian)
        inline int singularity_power(){ return 2; };

        // Use GKPY phase shift, smoothly extrapolated to pi 
        inline double phase_shift(double s){ return GKPY::phase_shift(1, 1, s); };

        // Since the prefactors are trivial, the kernel also is
        inline complex kernel(int iso_id, complex s, complex t)
        { 
            // Only P-wave allowed in the cross channel
            if (iso_id != V_to_3pi::kP_wave) return NaN<complex>();

            complex k  = _kinematics->kacser(s);
            complex kz = 2*t-_kinematics->Sigma()+s;
            return 3*(k*k - kz*kz); // We've multiplied by k^2 which is why singularity_power() = 2
        };

        // Choose default parameters for this isobar
        static const settings default_settings()
        {
            settings sets;
            sets._angular_integrator_depth    = 10;
            sets._dispersion_integrator_depth = 15;
            sets._interp_energy_low           = 2.;
            sets._interp_energy_high          = 20.;
            sets._interp_offset               = 0.2;
            sets._interp_points_low           = 100;
            sets._interp_points_high          = 50;

            sets._matching_intervals  = {0.03, 0.07, 0.08};
            sets._expansion_offsets   = {0.03, 0.1, 0.08};

            return sets;
        };
    };

}; /* namespace iterateKT */ }; // namespace V_to_3pi

#endif // V_TO_3PI_HPP