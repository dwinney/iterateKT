// Calculate partial widths and dalitz parameters for K -> 3π using parameters in [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2403.17570
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "plotter.hpp"

#include "isobars/kaon.hpp"
#include "amplitudes/kaon.hpp"
#include "K_3pi/data.hpp"

void calculate_observables()
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
                                 
    std::vector<complex> pars = {
        2.8147348,   // alpha_1
        -6418.4101,  // beta_1
        6080.5513,   // gamma_1
        -11550.85,   // zeta_1
        0.,          // eta
        -53.094403,  // alpha_3
        399.23876,   // beta_3
        -900.79741,  // gamma_3
        1071.8038,   // zeta_3
        -64.398136,  // mu
        4602.7463    // nu
    };
    amp->set_parameters(amp->process_fitter_parameters(pars));

    // --------------------------------------------------------------------------
    // When calculating widths integrate over the physical phase space, not isopin limit
    
    // Grab data
    auto P_ppm  = kaon::get_dalitz_data(option::P_ppm);
    auto P_zzp  = kaon::get_dalitz_data(option::P_zzp);
    auto L_pmz  = kaon::get_dalitz_data(option::L_pmz);
    auto L_zzz  = kaon::get_dalitz_data(option::L_zzz);
    auto lambda = kaon::get_lambda_data();
    std::vector<data_set> all_data = {P_ppm, P_zzp, L_pmz, L_zzz, lambda};
    
    double eps = 1E-3;
    std::array<double,5> dpars; // Dalitz plot parameters
    std::array<double,3> chi2_dpar;

    line(); divider();
    print<10,17>("", "Exp. value", "Our fit", "χ²");
    line();
    print<25>("", "K⁺ → π⁺π⁺π⁻"); divider(4);
    amp->set_option(option::P_ppm);
    dpars = amp->get_dalitz_parameters(eps);
    chi2_dpar = kaon::fit::chi2_dpars(P_ppm, amp);
    print<11,17>("Γ", "2.9590(218)",  physical_width(amp, option::P_ppm), kaon::fit::chi2_width(P_ppm, amp));
    print<10,17>("g", "-0.21134(17)", dpars[0], chi2_dpar[0]);
    print<10,17>("h", "0.0185(4)",    dpars[1], chi2_dpar[1]);
    print<10,17>("k", "-0.00463(14)", dpars[3], chi2_dpar[2]);
    line();

    print<25>("","K⁺ → π⁰π⁰π⁺"); divider(4);
    amp->set_option(option::P_zzp);
    dpars = amp->get_dalitz_parameters(eps);
    chi2_dpar = kaon::fit::chi2_dpars(P_zzp, amp);
    print<11,17>("Γ", "0.9438(150)", physical_width(amp, option::P_zzp), kaon::fit::chi2_width(P_zzp, amp));
    print<10,17>("g", "0.626(7)",    dpars[0], chi2_dpar[0]);
    print<10,17>("h", "0.052(8)",    dpars[1], chi2_dpar[1]);
    print<10,17>("k", "0.0054(35)",  dpars[3], chi2_dpar[2]);
    line();
    
    print<25>("","KL → π⁰π⁰π⁰"); divider(4);
    amp->set_option(option::L_zzz);
    dpars = amp->get_dalitz_parameters(eps);
    chi2_dpar = kaon::fit::chi2_dpars(L_zzz, amp);
    print<11,17>("Γ", "2.5417(352)", physical_width(amp, option::L_zzz), kaon::fit::chi2_width(L_zzz, amp));
    print<10,17>("h", "-0.0061(10)", dpars[1], chi2_dpar[1]);
    line();

    print<25>("","KL → π⁺π⁻π⁰"); divider(4);
    amp->set_option(option::L_pmz);
    dpars = amp->get_dalitz_parameters(eps);
    chi2_dpar = kaon::fit::chi2_dpars(L_pmz, amp);
    print<11,17>("Γ", "1.6200(102)", physical_width(amp, option::L_pmz), kaon::fit::chi2_width(L_pmz, amp));
    print<10,17>("g", "0.678(8)",    dpars[0], chi2_dpar[0]);
    print<10,17>("h", "0.076(6)",    dpars[1], chi2_dpar[1]);
    print<10,17>("k", "0.0099(15)",  dpars[3], chi2_dpar[2]);
    line();

    print<25>("","KS → π⁺π⁻π⁰"); divider(4);
    complex   lam = interference_lambda(amp);
    auto chi2_lam = kaon::fit::chi2_lambda(lambda, amp);
    amp->set_option(option::S_pmz);
    print<11,17>("Γ",    "0.0026(7)",   physical_width(amp, option::S_pmz), "    -");
    print<11,17>("Re λ", "0.0334(52)",  real(lam), chi2_lam[0]);
    print<11,17>("Im λ", "-0.0108(48)", imag(lam), chi2_lam[1]);
    line();

    print<8,19>("", "", "Total χ²:", kaon::fit::fcn(all_data, amp));

    divider();
};