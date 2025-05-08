// Fit KT amplitudes for K→3π decay data
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "plotter.hpp"
#include "fitter.hpp"

#include "isobars/kaon.hpp"
#include "amplitudes/kaon.hpp"
#include "data/kaon.hpp"

void fit()
{
    using namespace iterateKT;
    using iterateKT::complex;
    
    // --------------------------------------------------------------------------
    // Set up the amplitude from previously calculated isobars

    // Amplitude itself is given by the isospin limit
    kinematics kin = new_kinematics(M_KAON_AVG, M_PION_PM);
    amplitude  amp = new_amplitude<K_3Pi>(kin, "K -> 3π");

    // Empty array of subtraction indices for isobars with no polynomial
    std::vector<uint> empty = {};

    // Isobars for ΔI = 1/2 amplitude
    isobar M0 = amp->add_isobar<I1_S0>({0, 1, 2}, 2, id::dI1_I1_S0, "M0");
    isobar M1 = amp->add_isobar<I1_P1>({1},       1, id::dI1_I1_P1, "M1");
    isobar M2 = amp->add_isobar<I1_S2>(empty,     2, id::dI1_I1_S2, "M2"); 
    isobar G1 = amp->add_isobar<I0_P1>({1},       1, id::dI1_I0_P1, "G1");

    // and for the ΔI = 3/2 
    isobar N0 = amp->add_isobar<I1_S0>({0, 1, 2}, 2, id::dI3_I1_S0, "N0");
    isobar N1 = amp->add_isobar<I1_P1>({1},       1, id::dI3_I1_P1, "N1");
    isobar N2 = amp->add_isobar<I1_S2>(empty,     2, id::dI3_I1_S2, "N2"); 
    isobar H1 = amp->add_isobar<I2_P1>({0, 1},    1, id::dI3_I2_P1, "H1");
    isobar H2 = amp->add_isobar<I2_S2>(empty,     2, id::dI3_I2_S2, "H2");

    // Path to precalculated isobar files
    std::string path   = "/scripts/kaon/basis_functions/";
    std::string prefix = "K_3pi_";
    // Import everything 
    for (auto iso : amp->get_isobars()) iso->import_iteration<11>(path+prefix+iso->name()+".dat");

    // --------------------------------------------------------------------------
    // Set up fitter

    fitter<kaon::fit> fitter(amp);

    // Fit tends to be slow so its nice to have print level != 0 to see some progress
    fitter.set_print_level(4);
    fitter.set_tolerance(1);
    
    // Parameter labels and starting guess
    std::vector<std::string> labels = {"alpha_1", "beta_1", "gamma_1", "zeta_1", "eta",
                                       "alpha_3", "beta_3", "gamma_3", "zeta_3", "mu", "nu"};
    fitter.set_parameter_labels(labels);

    // Add data
    auto P_ppm = kaon::get_data(option::P_ppm);
    auto P_zzp = kaon::get_data(option::P_zzp);
    auto L_pmz = kaon::get_data(option::L_pmz);
    auto L_zzz = kaon::get_data(option::L_zzz);

    fitter.add_data(P_ppm);
    fitter.add_data(P_zzp);
    fitter.add_data(L_pmz);
    fitter.add_data(L_zzz);

    // We only fit real parts so force only fitting real parts
    for (auto par : labels) fitter.make_real(par);
    // Also fix eta which isnt fit
    fitter.fix_parameter("eta", 0.);
    
    // Initial guess from a previous fit
    std::vector<complex> initial = {
        58.338563,   // alpha_1
        -7070.9062,  // beta_1
        7011.0392,   // gamma_1
        -10312.59,   // zeta_1
        // 0.,       // eta (not fit)
        -47.263109,  // alpha_3
        323.41078,   // beta_3
        -730.62241,  // gamma_3
        1184.8828,   // zeta_3
        -74.541407,  // mu
        4674.9023    // nu
    };
    fitter.do_fit(initial);
};