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
#include "phaseshifts/GKPY.hpp"

namespace iterateKT { namespace V_to_3pi
{
    // Isobar id's
    static const int kP_wave = 0;
    
    // The P-wave is the dominant isobar
    class P_wave : public raw_isobar
    {
        // -----------------------------------------------------------------------
        public: 
        
        // Constructor 
        P_wave(kinematics xkin, int nsub, settings sets) : raw_isobar(xkin, nsub, sets)
        {
            _settings._angular_integrator_depth = 10;
        };

        // Virtual functions
        inline unsigned int  id() { return kP_wave;  };
        inline std::string name() { return "P-wave"; };

        // Because the P-wave involes a sintheta = 1-z^2, we have two power of 1/kappa
        // which lead to pseudo threshold singularities
        // The TOTAL singularity power is always +1 from this (one factor from the jacobian)
        inline int singularity_power(){ return 2; };

        // Use GKPY phase shift, smoothly extrapolated to pi 
        inline double phase_shift(double s){ return GKPY::phase_shift(1, 1, s); };

        // The prefactors are trivial for the J = 1
        inline complex prefactor(channel stu, complex s, complex t, complex u){ return 1.; };

        // Since the prefactors are trivial, the kernel also is
        inline complex kernel(int iso_id, complex s, complex t)
        { 
            // Only P-wave allowed in the cross channel
            if (iso_id != V_to_3pi::kP_wave) return NaN<complex>();

            complex k  = _kinematics->kacser(s);
            complex kz = 2*t-_kinematics->Sigma()+s;
            return 3*(k*k - kz*kz); // We've multiplied by k^2 which is why singularity_power() = 2
        };
    };

}; /* namespace iterateKT */ }; // namespace V_to_3pi

#endif // V_TO_3PI_HPP