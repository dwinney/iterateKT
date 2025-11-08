// Do a fit of K+ -> π+ π+ π- decay width and dalitz plot parameters
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "utilities.hpp"
#include "fitter.hpp"

#include "kaon_fitter.hpp"
#include "amplitudes/kaon.hpp"
#include "isobars/pseudoscalar.hpp"
#include "TRandom.h"

void fit()
{
    using namespace iterateKT;
    using iterateKT::complex;

    std::string method = "Simplex";
    double tolerance   = 1E-6;
    uint print_level   = 4;
    uint strategy      = 1;

    // Either choose N and leave initial_guess empty
    // (do N fits to do with randomized initial guesses)
    uint N             = 10;
    std::vector<complex> initial_guess = {};
    // or specify intitial_guess and do a single fit (ignore N)
    // std::vector<complex> initial_guess = {-0.726863775, -0.739814627, -4.89801606};

    // -----------------------------------------------------------------------
    
    // Set up general kinematics so everything knows masses
    // Assume masses are given in terms of pion mass
    kinematics kin = new_kinematics(M_KAON_PM, M_PION_PM);

    // Set up our amplitude 
    amplitude amp = new_amplitude<charged_kaon>(kin);
    amp->set_name("K⁺ → π⁺π⁺π⁻");

    settings sets = default_settings();
    sets._cauchy_integrator_depth  = 5;
    sets._pseudo_integrator_depth  = 5;
    
    // Add all the isobars, note the order they are added will be the order
    // the basis functions are generated
    std::vector<uint> empty = {}; // Pass empty to isobars with no sub polynomials
    isobar F0 = amp->add_isobar<I1_S0>(       2, id::I1_S0, "F0", sets); 
    isobar F1 = amp->add_isobar<I1_P1>(       1, id::I1_P1, "F1", sets); 
    isobar F2 = amp->add_isobar<I1_S2>(empty, 1, id::I1_S2, "F2", sets);
    isobar H1 = amp->add_isobar<I2_P1>(       1, id::I2_P1, "H1", sets); 
    isobar H2 = amp->add_isobar<I2_S2>(empty, 1, id::I2_S2, "H2", sets);

    // Path to precalculated isobar files
    std::string path   = "/scripts/kaon/basis_functions/basis_";
    // Import everything 
    for (auto iso : amp->get_isobars()) iso->import_iteration<6>(path+iso->name()+".dat");

    // -----------------------------------------------------------------------
    // Set up data
      
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
    
    // -----------------------------------------------------------------------
    // Set up fitter

    fitter<kaon::fit> fitter(amp, method, tolerance);
    fitter.set_print_level(print_level);
    fitter.set_strategy(strategy);

    // Add data from above
    fitter.add_data(dPtoPPM);
    fitter.add_data(dPtoZZP);

    // Set up parameters (all real)
    std::vector<std::string> labels = {"alpha", "beta", "gamma", "zeta"};
    fitter.set_parameter_labels(labels);
    for (auto par : labels) fitter.make_real(par);
    fitter.fix_parameter("zeta", 1);

    TRandom * guesser = new TRandom(0);
    double best_chi2 = -1; int status;
    std::vector<iterateKT::complex> best_pars;

    if (initial_guess.size() == 0)
    {
        for (int n = 0; n < N; n++)
        {
            initial_guess.clear();
            for (int i = 0; i < labels.size()-1; i++) initial_guess.push_back(guesser->Uniform(-1, 1));
            fitter.do_fit(initial_guess);
            if (best_chi2 == -1 || fitter.fcn() <= best_chi2)
            {
                status    = fitter.status();
                best_chi2 = fitter.fcn();
                best_pars = fitter.pars();
            };
        };
    }
    else
    {
        fitter.do_fit(initial_guess);
        status    = fitter.status();
        best_chi2 = fitter.fcn();
        best_pars = fitter.pars();
    };

    auto processed_pars = kaon::fit::process_fitter_parameters(best_pars, amp);
    amp->set_parameters(processed_pars);
    line();
    divider<20>(2);
    print("Status ", status);
    print("Best chi2 ", best_chi2);
    divider<20>(2);
    print<15,20>("alpha", processed_pars[0]);
    print<15,20>("beta",  processed_pars[1]);
    print<15,20>("gamma", processed_pars[2]);
    print<15,20>("zeta",  processed_pars[3]);
    divider<20>(2);
    line();
        
    amp->set_option(option::P_ppm);
    auto   cdpars = amp->get_dalitz_parameters(kaon::fit::derivative_h);
    divider<20>(4); 
    print<15,20>("", "Fit value", "Exp. value", "chi2");
    divider<20>(4); centered<20>(4, "K+ -> pi+ pi+ pi-"); divider<20>(4);
    print<15,20>("g",     cdpars[0], cg,   norm((cdpars[0]-cg)/cdg));
    print<15,20>("h",     cdpars[1], ch,   norm((cdpars[1]-ch)/cdh));
    print<15,20>("k",     cdpars[3], ck,   norm((cdpars[3]-ck)/cdk));
    
    amp->set_option(option::P_zzp);
    auto   ndpars = amp->get_dalitz_parameters(kaon::fit::derivative_h);
    divider<20>(4); centered<20>(4, "K+ -> pi0 pi0 pi+"); divider<20>(4);
    print<15,20>("g",     ndpars[0], ng,   norm((ndpars[0]-ng)/ndg));
    print<15,20>("h",     ndpars[1], nh,   norm((ndpars[1]-nh)/ndh));
    print<15,20>("k",     ndpars[3], nk,   norm((ndpars[3]-nk)/ndk));
    divider<20>(4); line();
};