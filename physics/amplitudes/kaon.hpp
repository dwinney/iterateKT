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
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "TMatrixD.h"

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
                case (id::dI1_I0_P1): return (s2-s3);
                case (id::dI3_I2_P1): return (s2-s3);
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
    
    enum class option : unsigned int { P_ppm, P_zzp, L_pmz, S_pmz, L_zzz };

    class K_3pi : public raw_amplitude
    {
        public:
        K_3pi(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {
            p_ppm = new_amplitude<Kp_PipPipPim>(xkin);
            p_zzp = new_amplitude<Kp_PizPizPip>(xkin);
            L_pmz = new_amplitude<KL_PipPimPiz>(xkin);
            S_pmz = new_amplitude<KS_PipPimPiz>(xkin);
            L_zzz = new_amplitude<KL_PizPizPiz>(xkin);
            current = p_ppm;
        };

        // We only fit the real parts of the parameters while the imaginary parts are
        // given by requiring Taylor invariants have vanishing imaginary parts
        inline std::vector<complex> process_fitter_parameters(std::vector<complex> in_pars)
        {
            double eps = 1E-5, s0 = _kinematics->s0();

            //------------------------------------------------------------------------
            // First we fix the imaginary parts of the M's and N's (total 3π I=1)
            isobar F0 = get_isobar(id::dI1_I1_S0);
            isobar F1 = get_isobar(id::dI1_I1_P1);
            isobar F2 = get_isobar(id::dI1_I1_S2);
            std::array<isobar,3> F = {F0, F1, F2};

            // Grab all the Taylor coefficients
            // First index is isospin, second is basis function ID
            std::array<std::array<complex,4>,3> A, B, C;
            for (uint i = 0; i <= 2; i++)
            {
                for (uint n = 0; n <= 3; n++)
                {
                    A[i][n] = F[i]->basis_function(n, 0);
                    B[i][n] = F[i]->basis_derivative<1>(n, 0, eps);
                    C[i][n] = F[i]->basis_derivative<2>(n, 0, eps)/2.;
                };
            };

            // Construct the 4x4 matrix of Taylor invariants
            // First index is which invariant, and second is basis_id
            TArrayD reT_data(16), imT_data(16);
            for (uint n = 0; n <= 3; n++)
            {
                // See first rows of Eq. 6.5 in [1]
                std::array<complex,4> T;
                T[0] = A[0][n] + s0*B[0][n] + 4*(A[2][n]+s0*B[2][n])/3;
                T[1] = 3*A[1][n] + B[0][n] - 5*B[2][n]/3 + 9*s0*(B[1][n] + 2*s0*C[1][n]);
                T[2] = 3*C[0][n] + 4*C[2][n];
                T[3] = C[2][n] + B[1][n] + 3*s0*C[1][n];
               
                for (int j = 0; j <= 3; j++)
                {
                    reT_data[4*j+n] = real(T[j]);  imT_data[4*j+n] = imag(T[j]);
                };
            };
            TMatrixD reT(4,4, reT_data.GetArray()), imT(4,4, imT_data.GetArray());

            // Now we actually solve the matrix equation relating reMu and imMu
            TMatrixD M = reT.Invert()*imT; M *= -1;

            Double_t reMu_data[4], reNu_data[4];
            for (int i = 0; i <= 3; i++)
            { 
                reMu_data[i] = real(in_pars[i]); reNu_data[i] = real(in_pars[i+5]);
            };
            TVectorD reMu(4, reMu_data), reNu(4, reNu_data);
            TVectorD imMu = M*reMu, imNu = M*reNu;
            
            //------------------------------------------------------------------------
            // Now do the same for the H's 

            // These inhabit basis functions 9 & 10
            std::array<isobar,3> H;
            H[0] = nullptr;
            H[1] = get_isobar(id::dI3_I2_P1); H[2] = get_isobar(id::dI3_I2_S2); 
            for (uint i = 1; i <= 2; i++)
            {
                for (uint n = 0; n <= 1; n++)
                {
                    A[i][n] = H[i]->basis_function(n+9, 0);
                    B[i][n] = H[i]->basis_derivative<1>(n+9, 0, eps);
                    C[i][n] = H[i]->basis_derivative<2>(n+9, 0, eps)/2.;
                };
            };

            // Construct the 4x4 matrix of Taylor invariants
            // First index is which invariant, and second is basis_id
            TArrayD reTp_data(4), imTp_data(4);
            for (uint n = 0; n <= 1; n++)
            {
                // See second rows of Eq. 6.5 in [1]
                std::array<complex,2> Tp;
                Tp[0] = 3*A[1][n] - B[2][n] + 9*s0*(B[1][n] + 2*C[1][n]);
                Tp[1] = 3*B[1][n] + C[2][n] + 9*s0*C[1][n];
                               
                for (int j = 0; j <= 1; j++)
                {
                    reTp_data[2*j+n] = real(Tp[j]);  imTp_data[2*j+n] = imag(Tp[j]);
                };
            };
            TMatrixD reTp(2,2, reTp_data.GetArray()), imTp(2,2, imTp_data.GetArray());
            TMatrixD Mp = reTp.Invert()*imTp; Mp *= -1;

            Double_t reNup_data[2];
            for (int i = 0; i <= 1; i++) reNup_data[i] = real(in_pars[i+9]);
            TVectorD reNup(2, reNup_data), imNup = Mp*reNup;

            //------------------------------------------------------------------------
            // Finally repackaged and return
            std::vector<complex> out_pars;
            for (int i = 0; i <= 3; i++) out_pars.push_back(reMu[i] +I*imMu[i]);
            out_pars.push_back(in_pars[4]); // Skip this one
            for (int i = 0; i <= 3; i++) out_pars.push_back(reNu[i] +I*imNu[i]);
            for (int i = 0; i <= 1; i++) out_pars.push_back(reNup[i]+I*imNup[i]);
            return out_pars;
        };

        inline double combinatorial_factor(){ return current->combinatorial_factor(); };
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
                case (option::P_ppm): current = p_ppm; return;
                case (option::P_zzp): current = p_zzp; return;
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

    // Because the kaon mass is so smal, the mass splittings between isospin projections
    // make a noticable difference in the phase-space.
    // So, even though we calculate KT isobars in isospin limit, we can integrate with the realistic 
    complex integrate_over_physical_phasespace(amplitude amp, std::array<option,2> opts, std::array<double,4> ms)
    {
        using namespace boost::math::quadrature;
        double mK = ms[0], m1 = ms[1], m2 = ms[2], m3 = ms[3];
        double sth = norm(m2+m3), pth = norm(mK-m1);
        double sigma = mK*mK + m1*m1 + m2*m2 + m3*m3;
        auto dGamma = [amp,sigma,mK,m1,m2,m3,opts](double s)
        {
            double kappa = sqrt(kallen(s, mK*mK, m1*m1))*sqrt(kallen(s,m2*m2,m3*m3))/s;
            double tmin  = (sigma - s - kappa)/2, tmax  = (sigma - s + kappa)/2;
            auto d2Gamma = [amp,s,opts](double t)
            {
                complex A, B;
                amp->set_option(opts[0]); A = amp->evaluate(s,t);
                if  (opts[0] == opts[1]) return conj(A)*A;
                amp->set_option(opts[1]); B = amp->evaluate(s,t);  
                return conj(A)*B;     
            };
            return gauss_kronrod<double,15>::integrate(d2Gamma, tmin, tmax, 0, 1.E-9, NULL);
        };
        return gauss_kronrod<double,15>::integrate(dGamma, sth, pth, 0, 1.E-9, NULL);
    };

    double physical_width(amplitude amp, option opt)
    {
        double mK, m1, m2, m3;
        switch (opt)
        {
            case (option::P_ppm): 
            {
                mK = M_KAON_PM; 
                m1 = M_PION_PM; m2 = M_PION_PM; m3 = M_PION_PM; 
                break;
            };
            case (option::L_zzz):
            {
                mK = M_KAON_AVG; 
                m1 = M_PION_0;  m2 = M_PION_0;  m3 = M_PION_0; 
                break;
            }; 
            case (option::P_zzp):
            {
                mK = M_KAON_PM; 
                m1 = M_PION_AVG; m2 = M_PION_AVG; m3 = M_PION_AVG; 
                break;
            };
            case (option::L_pmz):
            case (option::S_pmz):
            {
                mK = M_KAON_0; 
                m1 = M_PION_AVG; m2 = M_PION_AVG; m3 = M_PION_AVG; 
                break;
            };
            default: return NaN<double>();
        };

        double prefactors = 32*pow(2*PI*mK,3)*amp->combinatorial_factor();
        complex Gam = integrate_over_physical_phasespace(amp, {opt,opt}, {mK,m1,m2,m3});
        return real(Gam) / prefactors;
    };

    inline complex interference_lambda(amplitude amp)
    {
        double mK, m1, m2, m3;
        mK = M_KAON_0; 
        m1 = M_PION_AVG; m2 = M_PION_AVG; m3 = M_PION_AVG; 
        complex num = integrate_over_physical_phasespace(amp,{option::L_pmz,option::S_pmz}, {mK,m1,m2,m3});
        complex den = integrate_over_physical_phasespace(amp,{option::L_pmz,option::L_pmz}, {mK,m1,m2,m3});
        print(num, den);
        return num/den;
    };
};

#endif