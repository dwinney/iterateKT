// Implementation of the K -> π vector transition form factor
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Instituto de Ciencias Nucleares (ICN)
//               Universidad Autónoma Nacional de México (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] - https://arxiv.org/abs/2505.15309
// [2] - https://arxiv.org/abs/2510.17316
// ------------------------------------------------------------------------------

#ifndef KAON_PION_VECTOR_FORM_FACTOR_HPP
#define KAON_PION_VECTOR_FORM_FACTOR_HPP

#include "form_factor.hpp"

namespace iterateKT
{
    class K_to_pi_vector_TFF : public raw_form_factor
    {  
        public: 

        // Constructor
        K_to_pi_vector_TFF(uint nsub, amplitude amp, id proj, settings settings)
        : raw_form_factor(nsub, amp, proj, settings)
        {};

        // -----------------------------------------------------------------------
        // Virtual functions

        // The external current is given by the ππ vector form factor. 
        // We use the parameterization of [1] to incorporate ρ - ω mixing
        inline complex external_current(double s)
        {
            // Linear slope parameter 
            double alpha     = 0.0529;

            // ρ - ω mixing parameter
            double eps_omega = 0.156E-2;

            // ω pole given by the nominal PDG values
            double M_omega   = 782.66E-3;
            double G_omega   =   8.68E-3;

            complex P = 1 + alpha*s;
            complex G = 1 + eps_omega*s/(norm(M_omega)-s-I*G_omega*M_omega);

            // All the direct isobars assumed P-wave and thus all have the same omnes
            return P*G*_direct_isobars[0]->omnes(s);
        };

        // The normalizations are given in [2]
        inline double kinematic_factors(double s)
        {
            return PI * pow(_kinematics->rho(s), 3./2);
        };

        private:

        // 
    };
};

#endif