// Isobars relevant for the decay of meson with JP = 1-- into 3pi as in Ref. [1-2]
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] -  https://arxiv.org/abs/2006.01058
// [2] -  https://arxiv.org/abs/2304.09736
// ------------------------------------------------------------------------------

#ifndef VECTOR_ISOBARS_HPP
#define VECTOR_ISOBARS_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "GKPY.hpp"

namespace iterateKT
{ 
    // Ids for all our isobars. 
    // If we ignore rho-omega mixing, we have only one, else we have two
    enum class id : unsigned int { P_wave, charged, neutral, F_wave };

    // Choose default parameters for isobars
    inline settings default_settings()
    {
        settings sets;
        sets._infinitesimal           = 1E-8;
        sets._intermediate_energy     = 1.5;
        sets._cutoff                  = 20;
        sets._interpolation_offset    = 0.1;
        sets._interpolation_points    = {400, 10, 200};
        double xi_sth = 1E-3,  eps_sth = 1E-3;
        double xi_pth = 1E-3,  eps_pth = 1E-3;
        double xi_rth = 1E-3,  eps_rth = 1E-3;

        sets._exclusion_offsets   = {1E-1, 1E-1};
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};
        return sets;
    };

    // The P-wave is the dominant isobar
    // In terms of individual isobars this is the only one we need
    class P_wave : public raw_isobar
    {
        // -----------------------------------------------------------------------
        public: 
        
        // Constructor 
        P_wave(isobar_args args) : raw_isobar(args) {};
        
        // Because the P-wave involes a sintheta = 1-z^2, we have two power of 1/kappa
        // which lead to pseudo threshold singularities
        inline unsigned int angular_momentum(){ return 1; };

        // Use GKPY phase shift, smoothly extrapolated to pi 
        inline double phase_shift(double s){ return GKPY::phase_shift(1, 1, s); };

        // Kernels is 3*(1-z^2) but to remove kinematic singularities
        // we multiply by two powers of kappa
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            // Only P-wave allowed in the cross channel
            if (iso_id != id::P_wave) return 0.;
            complex k  = _kinematics->kacser(s), kz = _kinematics->kz(s,t);
            return 3*(k*k - kz*kz); // We've multiplied by k^2 which is why singularity_power() = 2
        };
    };

    class F_wave : public raw_isobar
    {
        // -----------------------------------------------------------------------
        public: 
        
        // Constructor 
        F_wave(isobar_args args) : raw_isobar(args) {};
        
        // Use phase of a simple BW for the rho3
        inline double phase_shift(double s)
        { 
            if (s <= _kinematics->sth()) return 0.;
            
            // Bare mass and width of the rho_3(1690)
            double m_rho3 = 1.688, gam_rho3 = 0.161;

            // Energy dependent width
            complex ps = _kinematics->momentum_final(s);
            complex pm = _kinematics->momentum_final(m_rho3*m_rho3);

            double r = 2; // GeV-1
            double z = norm(r*ps), z0 = norm(r*pm);
            double blatt_weisskopf = sqrt(z0*norm(z0-15) + 9*norm(2*z0-5))/sqrt(z*norm(z-15)+9*norm(2*z-5));

            double gamma =  gam_rho3*m_rho3/sqrt(s)*pow(real(ps/pm),7)*norm(blatt_weisskopf);

            complex BW = m_rho3*m_rho3/(m_rho3*m_rho3 - s - I*m_rho3*gamma);

            return arg(BW);
        };
        
        inline bool calculate_inhomogeneity(){ return false; };

        // We will ignore the inhomogeneity so singularity_power doesnt matter but it would be 6
        inline unsigned int angular_momentum(){ return 3; };
    };

    
    // ------------------------------------------------------------------------------
    // The following isobars are for the case that rho-omega mixing is handled explicitly

    class charged : public raw_isobar
    {
        // -----------------------------------------------------------------------
        public: 
        
        // Constructor 
        charged(isobar_args args) : raw_isobar(args) {};
        
        inline double phase_shift(double s){ return GKPY::phase_shift(1, 1, s); };
        inline unsigned int angular_momentum(){ return 1; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            complex k  = _kinematics->kacser(s), kz = _kinematics->kz(s,t);
            complex K11 = 3*(k*k - kz*kz);
            return (iso_id == id::charged || iso_id == id::neutral) ? K11/2 : 0.;
        };
    };

    class neutral : public raw_isobar
    {
        // -----------------------------------------------------------------------
        public: 
        
        // Constructor 
        neutral(isobar_args args) : raw_isobar(args) {};
        
        inline double phase_shift(double s){ return GKPY::phase_shift(1, 1, s); };
        inline unsigned int angular_momentum(){ return 1; };
        inline complex ksf_kernel(id iso_id, complex s, complex t)
        { 
            complex k  = _kinematics->kacser(s), kz = _kinematics->kz(s,t);
            complex K11 = 3*(k*k - kz*kz);
            return (iso_id == id::charged) ? K11 : 0;
        };
    };
}; // namespace iterateKT 

#endif // OMEGA_ISOBARS_HPP