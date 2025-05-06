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

        static double fcn(std::vector<data_set> & data_vector, amplitude to_fit)
        {
            double chi2 = 0;
            for (auto data : data_vector)
            {
                to_fit->set_option(data._option);
                double width = to_fit->width();
                auto   dpars = to_fit->get_dalitz_parameters(1.E-3, M_PION_PM*M_PION_PM);
                double g = dpars[0], h = dpars[1], k = dpars[3];

                // χ² from Γ
                chi2 += norm((width - data._z[0])/data._dz[0]); 
                // χ² from h
                chi2 += norm((h     - data._z[1])/data._dz[1]);

                if (data._type == kHOnly) continue;
                
                // χ² from g
                chi2 += norm((g     - data._z[2])/data._dz[2]);
                // χ² from k
                chi2 += norm((k     - data._z[3])/data._dz[3]);
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
                    0.0185,   // h
                    -0.21134, // g
                    -0.00463  // k
                };
                out._dz = {
                    218E-4,   // δΓ 
                    4E-4,     // δh
                    17E-5,    // δg
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
                    0.052,  // h
                    0.626,  // g
                    0.0054  // k
                };
                out._dz = {
                    150E-4, // δΓ 
                    8E-3,   // δh
                    7E-3,   // δg
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
                    0.076,  // h
                    0.678,  // g
                    0.0099  // k
                };
                out._dz = {
                    102E-4, // δΓ 
                    6E-3,   // δh
                    8E-3,   // δg
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