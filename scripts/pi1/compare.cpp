// Fit KT amplitudes for π1(1600) decay with only P-wave in [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2212.11767
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "plotter.hpp"
#include "fitter.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"
#include "COMPASS/data.hpp"

void compare()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // -----------------------------------------------------------------------
    // Operating options

    double m3pi = 1.4; // 3pi decay mass in GeV

    int Niter = 5; // Number of KT iterations

    // Best fit parameters for each of the three parameterizations
    complex omnes_norm              =  1379.5423*exp(I*-1.4768058);
    std::vector<complex> pars_1sub  = {1223.4514*exp(I*-1.4036273)};
    std::vector<complex> pars_2sub  = {34883.276*exp(I*-2.4690205),  59611.345*exp(I*0.69909863)};

    double smin = 0, smax = 1.5;

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution

    // Set up general kinematics so everything knows masses
    kinematics kin = new_kinematics(m3pi, M_PION);
    
    // Set up our amplitude 
    amplitude amp_1sub = new_amplitude<pi1>(kin);
    amp_1sub->add_isobar<P_wave>(1, id::P_wave);

    amplitude amp_2sub = new_amplitude<pi1>(kin);
    amp_2sub->add_isobar<P_wave>(2, id::P_wave);

    // Iterate
    divider();
    print("Solving KT with " + to_string(Niter) + " iterations:");
    line();
    timer timer; 
    timer.start();
    for (int i = 0; i < Niter; i++)
    {
        amp_1sub->iterate();
        amp_2sub->iterate();
        timer.lap("iteration " + to_string(i+1));
    };
    timer.stop();
    line();

    timer.print_elapsed();
    divider(); 

    // -----------------------------------------------------------------------
    // Plot results

    plotter plotter;

    plot p = plotter.new_plot();
    p.set_ranges({smin, smax}, {-5.5,10});
    p.set_legend(0.65, 0.75);
    p.set_labels("#sigma  [GeV^{2}]", "[ #it{f}#kern[-0.4]{_{1}}#kern[-1]{^{1}}(#sigma) - #it{f}#kern[-0.4]{_{1}}#kern[-1]{^{1}}(0) ] / 10^{3}");

    isobar pwave1 = amp_1sub->get_isobar(id::P_wave);
    isobar pwave2 = amp_2sub->get_isobar(id::P_wave);
    
    p.add_curve({smin, smax}, [&](double s){ return real(omnes_norm*pwave1->omnes(s+IEPS) - omnes_norm*pwave1->omnes(0))/1000; }, solid( jpacColor::Blue, "Omn#grave{e}s"));
    p.add_curve({smin, smax}, [&](double s){ return imag(omnes_norm*pwave1->omnes(s+IEPS) - omnes_norm*pwave1->omnes(0))/1000; }, dashed(jpacColor::Blue));

    amp_1sub->set_parameters(pars_1sub);
    p.add_curve({smin, smax}, [&](double s){ return real(pwave1->evaluate(s+IEPS) - pwave1->evaluate(0))/1000; }, solid( jpacColor::Red, "Unsubtracted"));
    p.add_curve({smin, smax}, [&](double s){ return imag(pwave1->evaluate(s+IEPS) - pwave1->evaluate(0))/1000; }, dashed(jpacColor::Red));
    
    amp_2sub->set_parameters(pars_2sub);
    p.add_curve({smin, smax}, [&](double s){ return real(pwave2->evaluate(s+IEPS) - pwave2->evaluate(0))/1000; }, solid( jpacColor::Green, "Once-subtracted"));
    p.add_curve({smin, smax}, [&](double s){ return imag(pwave2->evaluate(s+IEPS) - pwave2->evaluate(0))/1000; }, dashed(jpacColor::Green));

    p.save("compare.pdf");
   
};