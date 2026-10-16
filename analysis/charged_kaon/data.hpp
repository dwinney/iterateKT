// "Data sets" for the K → 3π decay to interface with fitters
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Instituto de Ciencias Nucleares (ICN)
//               Universidad Nacional Autonoma de Mexico (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#ifndef KAON_DATA_HPP
#define KAON_DATA_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "utilities.hpp"
#include "data_set.hpp"
#include "amplitudes/kaon.hpp"

namespace iterateKT { namespace kaon {

    // ID types, Dalitz parameterizations (Weinberg, Cabbibo-Isidori, or Empirical) and Width
    static const int kDalitz_W = 0, kDalitz_C = 1, kDalitz_E = 2, kWidth = 3;

    // Differentiate different experimental sources
    enum experiment
    {
        // These use the standard Weinberg expansion
        Charged_NA48_2007,  /* NA482:2007exu */
        Neutral_PDG_AVG,   
        Neutral_IHEP_2005,  /* Akopdzhanov:2005nb */
        Neutral_ISTRA_2002, /* Ajinenko:2002mg */
        // These use the Cabibbo-Isidori expansion
        Neutral_NA48_2009,  /* Batley:2009ubw */
        // This one uses an empirical parameterization
        Neutral_NA48_2010   /* NA482:2010gwp */
    };

    // The data_set class is really flexible and here we use it just as a method to 
    // pass all the dalitz plot parameters
    inline data_set dalitz_expansion(experiment exp)
    {
        option opt; 
        uint type;
        std::string id;
        double g, dg, h, dh, k, dk;         // Weinberg expansion parameters
        double a, da, b, db, p, dp, q, dq;  // extra parameters in Cabbibo-Isidori-like parameterization
        double s0;

        switch (exp)
        {
            /* NA482:2007exu */
            case experiment::Charged_NA48_2007: 
            { 
                opt  = option::P_ppm;
                type = kDalitz_W;
                id   = "Charged mode [NA48 2007]";
                g = -0.21134; dg = 1.7E-4;
                h = +0.01848; dh = 4.0E-4;
                k = -0.00463; dk = 1.4E-4; 
                s0 = (norm(M_KAON_PM) + 3*norm(M_PION_PM))/3;
                break;
            };
            
            case experiment::Neutral_PDG_AVG:
            {
                opt  = option::P_zzp;
                type = kDalitz_W;
                id   = "Neutral mode [PDG avg]";
                g = 0.626;  dg = 0.007;
                h = 0.052;  dh = 0.008;
                k = 0.0054; dk = 0.0035;
                s0 = (norm(M_KAON_PM) + norm(M_PION_PM) + 2*norm(M_PION_0))/3;    
                break; 
            };

            /* Akopdzhanov:2005nb */
            case experiment::Neutral_IHEP_2005:
            {
                opt  = option::P_zzp;
                type = kDalitz_W;
                id   = "Neutral mode [IHEP 2005]";
                g = 0.6259; dg = 4.3E-3 + 9.3E-3;         
                h = 0.0551; dh = 4.4E-3 + 8.6E-3;
                k = 0.0082; dk = 1.1E-3 + 1.4E-3;
                s0 = (norm(M_KAON_PM) + norm(M_PION_PM) + 2*norm(M_PION_0))/3;    
                break; 
            };

            /* Ajinenko:2002mg */
            case experiment::Neutral_ISTRA_2002:
            {
                opt  = option::P_zzp;
                type = kDalitz_W;
                id   = "Neutral mode [ISTRA 2002]";
                g = 0.627; dg = 0.004 + 0.010;
                h = 0.046; dh = 0.004 + 0.012;
                k = 0.001; dk = 0.001 + 0.002;
                s0 = (norm(M_KAON_PM) + norm(M_PION_PM) + 2*norm(M_PION_0))/3;    
                break; 
            };

            /* NA482:2010gwp */
            case experiment::Neutral_NA48_2010:
            {
                opt  = option::P_zzp;
                type = kDalitz_E;
                id   = "Neutral mode [NA48 2010]";
                g = 0.672;  dg = 0.011;
                h = -0.027; dh = 0.011;
                k = 0.0081; dk = 0.0005;
                a = -0.130; da = 0.007;
                b = -0.038; db = 0.009;
                p = 0.07;   dp = 0.03;
                q = 0.45;   dq = 0.06;
            };

            default: return error("kaon::get_dalitz_expansion(): Invalid id!", data_set());
        };

        data_set out;
        out._option = opt;
        out._type   = type;
        out._id     = id;
        // Instead of using the _x, _y, _z data, we just load everything into _extras
        out._extras["g"] = g; out._extras["dg"] = dg;
        out._extras["h"] = h; out._extras["dh"] = dh;
        out._extras["k"] = k; out._extras["dk"] = dk;
        out._extras["a"] = a; out._extras["da"] = da;
        out._extras["b"] = b; out._extras["db"] = db;
        out._extras["p"] = p; out._extras["dp"] = dp;
        out._extras["q"] = q; out._extras["dq"] = dq;
        // The physical s0 of the experimental process
        out._extras["s0"] = s0;
        
        return out;
    };

}; /* namespace charged_kaon */ }; /* namespace iterateKT */

#endif