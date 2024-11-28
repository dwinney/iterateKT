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

#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "GKPY.hpp"

namespace iterateKT
{ 
    // Ids for all our isobars, we only have one though
    enum class id : unsigned int { P_wave };

    // This defines the full amplitude, i.e. how the isobars are combined
    // Here is where we usually put the isospin combinations etc
    class vector : public raw_amplitude
    {
        public: 
        
        // Constructor
        vector(kinematics kin, std::string id) : raw_amplitude(kin,id)
        {};

        // We have no kinematic factors and only one isobar so simply return 1.
        // We're completely symmetric here so these are all the same
        inline complex prefactor_s(id iso_id, complex s, complex t, complex u){ return (iso_id == id::P_wave); };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u){ return prefactor_s(iso_id, t, s, u); };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u){ return prefactor_s(iso_id, u, t, s); };

        class P_wave;
    };

    
    // The P-wave is the dominant isobar
    // In terms of individual isobars this is the only one we need
    class vector::P_wave : public raw_isobar
    {
        // -----------------------------------------------------------------------
        public: 
        
        // Constructor 
        P_wave(kinematics xkin, subtractions subs, int maxsub, settings sets) 
        : raw_isobar(xkin, subs, maxsub, sets)
        {};

        // Virtual functions
        inline id get_id() { return id::P_wave;  };

        // Because the P-wave involes a sintheta = 1-z^2, we have two power of 1/kappa
        // which lead to pseudo threshold singularities
        // The TOTAL singularity power is always +1 from this (one factor from the jacobian)
        inline unsigned int singularity_power(){ return 2; };

        // Use GKPY phase shift, smoothly extrapolated to pi 
        inline double phase_shift(double s){ return GKPY::phase_shift(1, 1, M_PION*M_PION*s); };

        // Kernels is 3*(1-z^2) but to remove kinematic singularities
        // we multiply by two powers of kappa
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            // Only P-wave allowed in the cross channel
            if (iso_id != id::P_wave) return 0.;

            double m2  = _kinematics->m2(), M2 = _kinematics->M2();
            complex k  = _kinematics->kacser(s);
            complex kz = 2*t+s-M2-3*m2;
            return 3*(k*k - kz*kz); // We've multiplied by k^2 which is why singularity_power() = 2
        };

        // Choose default parameters for this isobar
        inline static const settings default_settings()
        {
            settings sets;
            sets._infinitesimal           = 1E-5;
            sets._intermediate_energy     = 60;
            sets._cutoff                  = 400;
            sets._interpolation_offset    = 1;
            sets._interpolation_points    = {400, 100, 200};
            double xi_sth = 0.3,  eps_sth = 0.5;
            double xi_pth = 0.3,  eps_pth = 0.3;
            double xi_rth = 0.8,  eps_rth = 1.5;
            sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
            sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};
            return sets;
        };
    };
}; // namespace iterateKT 

#endif // VECTOR_HPP