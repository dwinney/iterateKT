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
#include "K_3pi/data.hpp"

void fit()
{
    using namespace iterateKT;
    using iterateKT::complex;
    
    // --------------------------------------------------------------------------
    // Set up the amplitude from previously calculated isobars

    // Amplitude itself is given by the isospin limit
    kinematics kin = new_kinematics(M_KAON_AVG, M_PION_PM);
    amplitude  amp = new_amplitude<K_3pi>(kin, "K -> 3π");

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
    
    // Free parameters from Ref. [1] (ΔI = 1/2)
    complex alpha_1, beta_1, gamma_1, zeta_1, eta;     
    alpha_1 = +3.8;
    beta_1  = -676.1;
    gamma_1 = +559.7;
    zeta_1  = -1072.6;
    // ΔI = 3/2
    complex alpha_3, beta_3, gamma_3, zeta_3, mu, nu;  
    alpha_3 = -4.7;
    beta_3  = +26.7; 
    gamma_3 = -46.0;
    zeta_3  = +123.9;
    mu      = -2.04;
    nu      = +433.2;
    // Load parameters in the correct order (see order they were loaded above)
    std::vector<complex> initial = {alpha_1, beta_1, gamma_1, zeta_1,   
                                  alpha_3, beta_3, gamma_3, zeta_3, mu,  nu};
    fitter.do_fit(10*initial);
};