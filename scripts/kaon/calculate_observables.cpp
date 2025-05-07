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

void calculate_observables()
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

    // Free parameters from Ref. [1] (ΔI = 1/2)
    complex alpha_1, beta_1, gamma_1, zeta_1, eta;     
    alpha_1 = +3.8    - I*0.570;
    beta_1  = -676.1  + I*7.27;
    gamma_1 = +559.7  - I*16.80;
    zeta_1  = -1072.6 + I*7.57;
    eta     = +0.; // Undetermined assume its zero
    // ΔI = 3/2
    complex alpha_3, beta_3, gamma_3, zeta_3, mu, nu;  
    alpha_3 = -4.7    - I*2.37E-2;
    beta_3  = +26.7   + I*0.30; 
    gamma_3 = -46.0   - I*0.74;
    zeta_3  = +123.9  - I*0.28;
    mu      = -2.04   + I*4.8E-4;
    nu      = +433.2  - I*6.0E-4;

    // Load parameters in the correct order (see order they were loaded above)
    std::vector<complex> pars = {alpha_1, beta_1, gamma_1, zeta_1, eta,   
                                 alpha_3, beta_3, gamma_3, zeta_3, mu,  nu};
    amp->set_parameters(10*pars);
    
    // --------------------------------------------------------------------------
    // When calculating widths integrate over the physical phase space, not isopin limit
    
    double eps = 1E-3;
    std::array<double,5> dpars; // Dalitz plot parameters

    line(); divider();
    print("Observable", "From [1]", "Ours");
    line();
    print("","K⁺ → π⁺π⁺π⁻"); divider(3);
    amp->set_option(option::P_ppm);
    dpars = amp->get_dalitz_parameters(eps);
    
    print<15>("Width",  2.9865,  width_with_physical_masses(amp, option::P_ppm));
    print<15>("g", -0.21134, dpars[0]);
    print<15>("h",  0.0185,  dpars[1]);
    print<15>("k", -0.00464, dpars[3]);
    line();

    print("","K⁺ → π⁰π⁰π⁺"); divider(3);
    amp->set_option(option::P_zzp);
    dpars = amp->get_dalitz_parameters(eps);
    print<15>("Width", 0.8984,  width_with_physical_masses(amp, option::P_zzp));
    print<15>("g", 0.625, dpars[0]);
    print<15>("h", 0.058, dpars[1]);
    print<15>("k", 0.011, dpars[3]);
    line();
    
    print("","KL → π⁺π⁻π⁰"); divider(3);
    amp->set_option(option::L_pmz);
    dpars = amp->get_dalitz_parameters(eps);
    print<15>("Width", 1.6185, width_with_physical_masses(amp, option::L_pmz));
    print<15>("g", 0.675, dpars[0]);
    print<15>("h", 0.082, dpars[1]);
    print<15>("k", 0.011, dpars[3]);
    line();

    print("","KL → π⁰π⁰π⁰"); divider(3);
    amp->set_option(option::L_zzz);
    dpars = amp->get_dalitz_parameters(eps);
    print<15>("Width", 2.5532,  width_with_physical_masses(amp, option::L_zzz));
    print<15>("h", -0.00644, dpars[1]);
    line(); divider();
};