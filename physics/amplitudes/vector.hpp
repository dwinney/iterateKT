// Amplitudes relevant for the decay of meson with JP = 1-- into 3pi as in Ref. [1]
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

#ifndef VECTOR_AMPLITUDES_HPP
#define VECTOR_AMPLITUDES_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "GKPY.hpp"

#include"isobars/vector.hpp"

namespace iterateKT
{ 
    // This defines the full amplitude, i.e. how the isobars are combined
    // Here is where we usually put the isospin combinations etc
    class vector_decay : public raw_amplitude
    {
        public: 
        
        // Constructor
        vector_decay(kinematics kin, std::string id) : raw_amplitude(kin,id)
        {};

        // Divide by 3 polarizations
        inline double combinatorial_factor(){ return 3; };
        
        // We have no kinematic factors and only one isobar so simply return 1.
        // We're completely symmetric here so these are all the same
        inline complex prefactor_s(id iso_id, complex s, complex t, complex u){ return (iso_id == id::P_wave) ? csqrt(_kinematics->kibble(s,t,u))/2 : 0; };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u){ return prefactor_s(iso_id, t, s, u); };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u){ return prefactor_s(iso_id, u, t, s); };
    };

    // This amplitude allows neutral and charged channels to be implemented independently
    class rho_omega_mixing : public raw_amplitude
    {
        public: 
        
        // Constructor
        rho_omega_mixing(kinematics kin, std::string id) : raw_amplitude(kin,id)
        {};

        inline double helicity_factor(){ return 3; };
        
        // Kinematic factors are still symmetric
        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        { 
            return (iso_id == id::neutral || iso_id == id::charged) ? csqrt(_kinematics->kibble(s,t,u))/2 : 0; 
        };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u){ return prefactor_s(iso_id, t, s, u); };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u){ return prefactor_s(iso_id, u, t, s); };
    };
}; // namespace iterateKT 

#endif // OMEGA_AMPLITUDES_HPP