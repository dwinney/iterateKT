// Methods to interface K → 3π decay data with fitters
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef KAON_DATA_HPP
#define KAON_DATA_HPP

#include "constants.hpp"
#include "amplitudes/kaon.hpp"
#include "data_set.hpp"

namespace iterateKT { namespace kaon
{
    // Specify the fitter interface
    struct fit
    {
        // Static identifiers for data_set types
        // Have two 'types' one which contains only Gamma and h
        // and another with Gamma, g, h, k
        static const uint kAll = 0, kHOnly = 1;
        
        // String letting us know what is being fit
        static std::string data_type(int i)
        {
            switch (i)
            {
                case kAll:   return "Γ & {g, h, k}";
                case kHOnly: return "Γ & h";
                default: return "ERROR!";
            };
        };

        // Function being minimized, sum of chi2s of individual data sets and observables
        static double fcn(std::vector<data_set> & data_vector, amplitude to_fit)
        {
            double chi2_tot = 0;
            for (auto data : data_vector) chi2_tot += chi2(data, to_fit);
            return chi2_tot;
        };

        // Indidivudal chi2 from a single data set
        static double chi2(const data_set & data, amplitude to_fit)
        {
            to_fit->set_option(data._option);
            
            // output
            double chi2 = 0;

            // χ² from Γ
            chi2 += chi2_width(data, to_fit); 

            // χ² from g h k
            auto chi2_ghk = chi2_dpars(data, to_fit);
            for (auto chi2_i : chi2_ghk) chi2 += chi2_i;
            return chi2;
        };

        // Compare widths
        static double chi2_width(const data_set & data, amplitude to_fit)
        {
            double gam_th = width_with_physical_masses(to_fit, data._option);
            double gam_ex = data._z[0], dgam_ex = data._dz[0];
            return norm((gam_th - gam_ex)/dgam_ex);
        };

        // Compare g, h, k
        static std::array<double,3> chi2_dpars(const data_set & data, amplitude to_fit)
        {
            auto dpars = to_fit->get_dalitz_parameters(1.E-3, M_PION_PM*M_PION_PM);
            double g_th = dpars[0], h_th = dpars[1], k_th = dpars[3];
     
            std::array<double,3> chi2;
            if (data._type == kHOnly)
            {
                double  h_ex = data._z[1],  dh_ex = data._dz[1];
                chi2[1] = norm((h_th - h_ex)/dh_ex);
                chi2[0] = 0; chi2[2] = 0;
            }
            else
            {
                double  g_ex = data._z[1],   h_ex = data._z[2],   k_ex = data._z[3];
                double dg_ex = data._dz[1], dh_ex = data._dz[2], dk_ex = data._dz[3];
                chi2[0] = norm((g_th - g_ex)/dg_ex);
                chi2[1] = norm((h_th - h_ex)/dh_ex);
                chi2[2] = norm((k_th - k_ex)/dk_ex);
            };
            return chi2;
        };
    };

    inline data_set get_data(option opt)
    {
        data_set out;
        switch (opt)
        {
            case option::P_ppm:
            {
                out._id = "K+ -> pi+ pi+ pi-";
                out._z  = {
                    2.9590,   // Γ
                    -0.21134, // g
                    0.0185,   // h
                    -0.00463  // k
                };
                out._dz = {
                    218E-4,   // δΓ 
                    17E-5,    // δg
                    4E-4,     // δh
                    14E-5     // δk
                };
                out._option = option::P_ppm;
                out._type = fit::kAll;
                out._N    = 4;
                break;
            };
            case option::P_zzp:
            {
                out._id = "K+ -> pi0 pi0 pi+";
                out._z  = {
                    0.9438, // Γ
                    0.626,  // g
                    0.052,  // h
                    0.0054  // k
                };
                out._dz = {
                    150E-4, // δΓ 
                    7E-3,   // δg
                    8E-3,   // δh
                    35E-4   // δk
                };
                out._option = option::P_zzp;
                out._type = fit::kAll;
                out._N    = 4;
                break;
            };
            case option::L_pmz:
            {
                out._id = "KL -> pi+ pi- pi0";
                out._z  = {
                    1.6200, // Γ
                    0.678,  // g
                    0.076,  // h
                    0.0099  // k
                };
                out._dz = {
                    102E-4, // δΓ 
                    8E-3,   // δg
                    6E-3,   // δh
                    15E-4   // δk
                };
                out._option = option::L_pmz;
                out._type = fit::kAll;
                out._N    = 4;
                break;
            };
            case option::L_zzz:
            {
                out._id = "KL -> pi0 pi0 pi0";
                out._z = {
                    2.5417,  // Γ
                    -0.0061  // h
                };
                out._dz = {
                    352E-4,  // δΓ
                    10E-4    // δh
                };
                out._option = option::L_zzz;
                out._type = fit::kHOnly;
                out._N    = 2;
                break;
            };
            default: return out;
        };
        return out;
    };
}; /* namespace iterateKT */ }; /* namespace kaon_decay */

#endif