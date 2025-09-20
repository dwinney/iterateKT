// Calculate the imaginary parts of subtraction coefficients given
// their real parts assuming Taylor invariants are purely real
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
#include "data_set.hpp"

#include "isobars/pseudoscalar.hpp"
#include "amplitudes/kaon.hpp"
#include "kaon_fitter.hpp"

void compare()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::option;

    // -----------------------------------------------------------------------
    // Set up data
    
    double cgam, cdgam;
    cgam =  2.9590,  cdgam = 218E-4;
    data_set wPtoPPM;
    wPtoPPM._id     = "K+ -> pi+ pi+ pi-";
    wPtoPPM._z      = {cgam};
    wPtoPPM._dz     = {cgam};
    wPtoPPM._N      = 1;
    wPtoPPM._option = option::P_ppm;
    wPtoPPM._type   = kaon::fit::kWidth;
    
    // K+ -> pi+ pi+ pi- dalitz parameters
    double cg, cdg, ch, cdh, ck, cdk;
    cg   = -0.21134, cdg   = 17E-5;
    ch   =  0.0185,  cdh   = 40E-5;
    ck   = -0.00463, cdk   = 14E-5;

    data_set dPtoPPM;
    dPtoPPM._id     = "K+ -> pi+ pi+ pi-";
    dPtoPPM._z      = {cg,  ch,  ck};
    dPtoPPM._dz     = {cdg, cdh, cdk};
    dPtoPPM._N      = 3;
    dPtoPPM._option = option::P_ppm;
    dPtoPPM._type   = kaon::fit::kDalitz;

    // K+ -> pi0 pi0 pi+ dalitz parameters
    double ngam, ndgam, ng, ndg, nh, ndh, nk, ndk;
    ng   =  0.626,  ndg   = 70E-4;
    nh   =  0.052,  ndh   = 80E-4;
    nk   =  0.0054, ndk   = 35E-4;

    data_set dPtoZZP;
    dPtoZZP._id     = "K+ -> pi0 pi0 pi+";
    dPtoZZP._z      = { ng,  nh,  nk};
    dPtoZZP._dz     = {ndg, ndh, ndk};
    dPtoZZP._N      = 3;
    dPtoZZP._option = option::P_zzp;
    dPtoZZP._type   = kaon::fit::kDalitz;

    // --------------------------------------------------------------------------
    // Set up the amplitude from previously calculated isobars

    // Amplitude itself is given by the isospin limit
    kinematics kin = new_kinematics(M_KAON_AVG, M_PION_PM);
    amplitude  amp = new_amplitude<charged_kaon>(kin, "K -> 3π");

    // Empty array of subtraction indices for isobars with no polynomial
    std::vector<uint> empty = {};

    // Isobars for ΔI = 1/2 amplitude
    isobar F0 = amp->add_isobar<I1_S0>({0, 1, 2}, 2, id::I1_S0, "F0");
    isobar F1 = amp->add_isobar<I1_P1>({1},       1, id::I1_P1, "F1");
    isobar F2 = amp->add_isobar<I1_S2>(empty,     2, id::I1_S2, "F2"); 
    isobar H1 = amp->add_isobar<I2_P1>({0, 1},    1, id::I2_P1, "H1"); 
    isobar H2 = amp->add_isobar<I2_S2>(empty,     2, id::I2_S2, "H2");

    // Path to precalculated isobar files
    std::string path   = "/scripts/kaon/basis_functions/basis_";
    for (auto iso : amp->get_isobars()) iso->import_iteration<6>(path+iso->name()+".dat");

    // Free parameters from Ref. [1] (ΔI = 1/2)
    complex alpha_1, beta_1, gamma_1, zeta_1, eta;     
    alpha_1 = +3.8    - I*0.570;
    beta_1  = -676.1  + I*7.27;
    gamma_1 = +559.7  - I*16.80;
    zeta_1  = -1072.6 + I*7.57;
    // ΔI = 3/2
    complex alpha_3, beta_3, gamma_3, zeta_3, mu, nu;  
    alpha_3 = -4.7    - I*2.37E-2;
    beta_3  = +26.7   + I*0.30; 
    gamma_3 = -46.0   - I*0.74;
    zeta_3  = +123.9  - I*0.28;

    // I_3pi = 2
    mu      = -2.04   + I*4.8E-4;
    nu      = +433.2  - I*6.0E-4;

    // Load parameters in the correct order (see order they were loaded above)
    std::vector<complex> pars1 = {alpha_1+alpha_3, beta_1+beta_3, gamma_1+gamma_3, zeta_1+zeta_3, mu, nu};  
    std::vector<complex> pars2 = (kaon::fit::process_fitter_parameters(pars1, amp));
    std::vector<std::string> labels = {"alpha", "beta", "gamma", "zeta", "mu", "nu"};

    line(); divider();
    print<15>("par", "Re(par)", "Im(par) [1]", "Im(par) [Us]");
    divider(4);
    for (int i = 0; i < pars1.size(); i++) print<15>(labels[i], real(pars1[i]), imag(pars1[i]), imag(pars2[i]));
    divider(); line();

    amp->set_parameters(pars2);
    amp->set_option(option::P_ppm);
    auto   cdpars = amp->get_dalitz_parameters(kaon::fit::derivative_h);
    divider<20>(4); 
    print<15,20>("", "Fit value", "Exp. value", "chi2");
    divider<20>(4); centered<20>(4, "K+ -> pi+ pi+ pi-"); divider<20>(4);
    print<15,20>("g",     cdpars[0], cg,   norm((cdpars[0]-cg)/cdg));
    print<15,20>("h",     cdpars[1], ch,   norm((cdpars[1]-ch)/cdh));
    print<15,20>("k",     cdpars[3], ck,   norm((cdpars[3]-ck)/cdk));
    print<15,20>("j",     cdpars[2], 0);
    
    amp->set_option(option::P_zzp);
    auto   ndpars = amp->get_dalitz_parameters(kaon::fit::derivative_h);
    divider<20>(4); centered<20>(4, "K+ -> pi0 pi0 pi+"); divider<20>(4);
    print<15,20>("g",     ndpars[0], ng,   norm((ndpars[0]-ng)/ndg));
    print<15,20>("h",     ndpars[1], nh,   norm((ndpars[1]-nh)/ndh));
    print<15,20>("k",     ndpars[3], nk,   norm((ndpars[3]-nk)/ndk));
    print<15,20>("j",     ndpars[2], 0);
    divider<20>(4); line();
};
