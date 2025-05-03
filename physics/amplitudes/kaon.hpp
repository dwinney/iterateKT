// Isobars relevant for the decay of isospin-1/2 decay to 3pi in [1]
//
// One may notice that these isobars follow the same structure as those of 
// the eta -> 3pi (e.g. in "eta.hpp")
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

#ifndef KAON_AMPLITUDES_HPP
#define KAON_AMPLITUDES_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shift.hpp"
#include "isobars/kaon.hpp"

#include "Math/IntegratorMultiDim.h"
#include "Math/IntegrationTypes.h"

// The amplitudes are named with respect to isospin projections of the decay particle
// into three pions. The order matters in that the definitions of s, t, and u.

// For a general K_Pi1Pi2Pi3, we have:
// s = (K - Pi1)^2 = (Pi2 + Pi3)^2 = s1
// t = (K - Pi2)^2 = (Pi1 + Pi3)^2 = s2
// u = (K - Pi3)^2 = (Pi1 + Pi2)^2 = s3

namespace iterateKT
{
    // Different masses for isospin projectionss
    const double M_KAON_PM  = 0.493677;
    const double M_KAON_0   = 0.497611;
    const double M_KAON_AVG = (M_KAON_PM + M_KAON_0)/2;
    const double M_PION_PM  = 0.13957039;
    const double M_PION_0   = 0.1349768;
    const double M_PION_AVG = (M_PION_PM + M_PION_0)/2;

    //--------------------------------------------------------------------------
    // K⁺ → π⁺ π⁺ π⁻         
    class Kp_PipPipPim : public raw_amplitude
    {
        public: 
        Kp_PipPipPim(kinematics xkin, std::string id) : raw_amplitude(xkin, id){};
        inline double combinatorial_factor(){ return 2.; }; // 2 identical particles 
        inline complex prefactor_s(id iso_id, complex s1, complex s2, complex s3)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_S0): return -1;
                case (id::dI1_I1_P1): return +(s3-s2);
                case (id::dI1_I1_S2): return -1./3;
                case (id::dI3_I1_S0): return -1;
                case (id::dI3_I1_P1): return +(s3-s2);
                case (id::dI3_I1_S2): return -1./3;
                case (id::dI3_I2_P1): return +3*(s3-s2)/2;
                case (id::dI3_I2_S2): return +1./2;
                default: return 0;
            };
        };
        // factors with s2 (same as s1 factors with s1 <-> s2)
        inline complex prefactor_t(id iso_id, complex s1, complex s2, complex s3)
        { return prefactor_s(iso_id, s2, s1, s3); };
        // factors with s3
        inline complex prefactor_u(id iso_id, complex s1, complex s2, complex s3)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_S2): return -2;
                case (id::dI3_I1_S2): return -2;
                case (id::dI3_I2_S2): return -1;
                default: return 0;
            };
        };
    };

    //--------------------------------------------------------------------------
    // K+ -> pi0 pi0 pi+         
    class Kp_PizPizPip : public raw_amplitude
    {
        public: 
        Kp_PizPizPip(kinematics xkin, std::string id) : raw_amplitude(xkin, id){};
        inline double combinatorial_factor(){ return 2.; }; // 2 identical particles
        inline complex prefactor_s(id iso_id, complex s1, complex s2, complex s3)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_P1): return +(s3-s2);
                case (id::dI1_I1_S2): return +1;
                case (id::dI3_I1_P1): return +(s3-s2);
                case (id::dI3_I1_S2): return +1;
                case (id::dI3_I2_P1): return -3*(s3-s2)/2;
                case (id::dI3_I2_S2): return -1./2;
                default: return 0;
            };
        };
        inline complex prefactor_t(id iso_id, complex s1, complex s2, complex s3)
        { return prefactor_s(iso_id, s2, s1, s3); };

        inline complex prefactor_u(id iso_id, complex s1, complex s2, complex s3)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_S0): return +1;
                case (id::dI1_I1_S2): return -2./3;
                case (id::dI3_I1_S0): return +1;
                case (id::dI3_I1_S2): return -2./3;
                case (id::dI3_I2_S2): return +1;
                default: return 0;
            };
        };
    };

    
    //--------------------------------------------------------------------------
    // KL -> pi+ pi- pi0         
    class KL_PipPimPiz : public raw_amplitude
    {
        public: 
        KL_PipPimPiz(kinematics xkin, std::string id) : raw_amplitude(xkin, id){};
        inline complex prefactor_s(id iso_id, complex s1, complex s2, complex s3)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_P1): return -(s3-s2);
                case (id::dI1_I1_S2): return -1;
                case (id::dI3_I1_P1): return +2*(s3-s2);
                case (id::dI3_I1_S2): return +2;

                default: return 0;
            };
        };
        inline complex prefactor_t(id iso_id, complex s1, complex s2, complex s3)
        { return prefactor_s(iso_id, s2, s1, s3); };
        inline complex prefactor_u(id iso_id, complex s1, complex s2, complex s3)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_S0): return -1;
                case (id::dI1_I1_S2): return +2./3;
                case (id::dI3_I1_S0): return +2;
                case (id::dI3_I1_S2): return -4./3;
                default: return 0;
            };
        };
    };

    //--------------------------------------------------------------------------
    // KS -> pi+ pi- pi0         
    class KS_PipPimPiz : public raw_amplitude
    {
        public: 
        KS_PipPimPiz(kinematics xkin, std::string id) : raw_amplitude(xkin, id){};
        inline complex prefactor_s(id iso_id, complex s1, complex s2, complex s3)
        {
            switch (iso_id)
            {
                case (id::dI1_I0_P1): return -(s3-s2);
                case (id::dI3_I2_P1): return -(s3-s2);
                case (id::dI3_I2_S2): return +1;
                default: return 0;
            };
        };
        inline complex prefactor_t(id iso_id, complex s1, complex s2, complex s3)
        { return - prefactor_s(iso_id, s2, s1, s3); }; // Get a minus sign
        inline complex prefactor_u(id iso_id, complex s1, complex s2, complex s3)
        {
            switch (iso_id)
            {
                case (id::dI1_I0_P1): return +(s1-s2);
                case (id::dI3_I2_P1): return -2*(s1-s2);
                default: return 0;
            };
        };
    };

    
    //--------------------------------------------------------------------------
    // KL -> pi0 pi0 pi0         
    class KL_PizPizPiz : public raw_amplitude
    {
        public: 
        KL_PizPizPiz(kinematics xkin, std::string id) : raw_amplitude(xkin, id){};     
        inline double combinatorial_factor(){ return 6.; }; // 3 identical particles
        inline complex prefactor_s(id iso_id, complex s1, complex s2, complex s3)
        {
            switch (iso_id)
            {
                case (id::dI1_I1_S0): return +1;
                case (id::dI1_I1_S2): return +4./3;
                case (id::dI3_I1_S0): return -2;
                case (id::dI3_I1_S2): return -8./3;
                default: return 0;
            };
        };
        inline complex prefactor_t(id iso_id, complex s1, complex s2, complex s3)
        { return prefactor_s(iso_id, s2, s1, s3); };
        inline complex prefactor_u(id iso_id, complex s1, complex s2, complex s3)
        { return prefactor_s(iso_id, s3, s2, s1); };
    };

    //--------------------------------------------------------------------------
    // Generic K -> 3 pi this will contain all the above which can be accessed via the set_option() function
    
    enum class option : unsigned int { p_ppm, p_zzp, L_pmz, S_pmz, L_zzz };

    class K_3Pi : public raw_amplitude
    {
        public:
        K_3Pi(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {
            p_ppm = new_amplitude<Kp_PipPipPim>(xkin);
            p_zzp = new_amplitude<Kp_PizPizPip>(xkin);
            L_pmz = new_amplitude<KL_PipPimPiz>(xkin);
            S_pmz = new_amplitude<KS_PipPimPiz>(xkin);
            L_zzz = new_amplitude<KL_PizPizPiz>(xkin);
            current = p_ppm;
        };
    
        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        { return current->prefactor_s(iso_id, s, t, u); };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        { return current->prefactor_t(iso_id, s, t, u); };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        { return current->prefactor_u(iso_id, s, t, u); };

        inline void set_option(option opt)
        {
            switch (opt)
            {
                case (option::p_ppm): current = p_ppm; return;
                case (option::p_zzp): current = p_zzp; return;
                case (option::L_pmz): current = L_pmz; return;
                case (option::S_pmz): current = S_pmz; return;
                case (option::L_zzz): current = L_zzz; return;
                default: return;
            };
        };

        private:
        amplitude current; 
        amplitude p_ppm, p_zzp, L_pmz, S_pmz, L_zzz;
    };

    // Because the kaon mass is so small, the mass splittings between isospin projections
    // make a noticable difference in the phase-space.
    // So, even though we calculate KT isobars in isospin limit, we can integrate with the realistic 
    double width_with_physical_masses(amplitude amp, double mK, std::array<double,3> m)
    {
        double m1 = m[0], m2 = m[1], m3 = m[2];
        double sth = norm(m2+m3), pth = norm(mK-m1);
        double prefactors = 32*pow(2*PI*mK,3)*amp->combinatorial_factor();

        // Doubly differential width with given masses not those in ampitude::_kinematics
        auto d2Gamma = [&] (const double * sz)
        {
            double s = sz[0], z = sz[1];
            double sigma = mK*mK + m1*m1 + m2*m2 + m3*m3;
            double kappa = sqrt(kallen(s, mK*mK, m1*m1))*sqrt(kallen(s,m2*m2,m3*m3))/s;
            double t     = (sigma - s - kappa*z)/2;
            return kappa/2 * norm(amp->evaluate(s,t)) / prefactors;
        };

        ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::Type::kADAPTIVE, 1E-5, 1E-5);
        double min[2] = {sth, -1.}, max[2] = {pth, 1.};
        return ig.Integral(d2Gamma, 2, min, max);
    };
};

#endif