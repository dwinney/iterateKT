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
        P_wave(kinematics xkin, int nsub) : raw_isobar(xkin, nsub)
        {};

        // Virtual functions
        inline unsigned int  id() { return kP_wave;  };
        inline std::string name() { return "P-wave"; };

        // Use GKPY phase shift, smoothly extrapolated to pi 
        inline double phase_shift(double s){ return GKPY::phase_shift(1, 1, s); };

        // The prefactors are trivial for the J = 1
        inline complex prefactor(channel stu, complex s, complex t, complex u){ return 1.; };
    };

}; /* namespace iterateKT */ }; // namespace V_to_3pi

#endif // V_TO_3PI_HPP