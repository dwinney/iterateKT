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
    complex omnes_norm              =  630.32738;
    std::vector<complex> pars_1sub  = {563.1772};
    std::vector<complex> pars_2sub  = {10954.435,  18585.493*exp(I*3.1524194)};
    
    double smin = 0, smax = 2.5;

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution

    // Set up general kinematics so everything knows masses
    kinematics kin = new_kinematics(m3pi, M_PION);
    
    // Set up our amplitude 
    amplitude amp_1sub = new_amplitude<pi1>(kin);
    isobar pwave1 = amp_1sub->add_isobar<P_wave>(1, id::P_wave);

    amplitude amp_2sub = new_amplitude<pi1>(kin);
    isobar pwave2 = amp_2sub->add_isobar<P_wave>(2, id::P_wave);

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

    // Set parameters
    amp_1sub->set_parameters(pars_1sub);
    amp_2sub->set_parameters(pars_2sub);

    // -----------------------------------------------------------------------
    // Plot results

    plotter plotter;

    plot p1 = plotter.new_plot();
    p1.set_curve_points(400);
    p1.set_ranges({smin, smax}, {-3.5,5});
    p1.set_legend(0.68, 0.75);
    p1.set_labels("#sigma  [GeV^{2}]", "( #it{f}(#sigma) - #it{f}(0) ) / 10^{3}");

    p1.add_curve({smin, smax}, [&](double s){ return real(omnes_norm*pwave1->omnes(s+IEPS) - omnes_norm*pwave1->omnes(0))/1E3; }, solid( jpacColor::Blue, "Omnes"));
    p1.add_curve({smin, smax}, [&](double s){ return imag(omnes_norm*pwave1->omnes(s+IEPS) - omnes_norm*pwave1->omnes(0))/1E3; }, dashed(jpacColor::Blue));
    p1.add_curve({smin, smax}, [&](double s){ return real(pwave1->evaluate(s+IEPS) - pwave1->evaluate(0))/1E3; },                 solid( jpacColor::Red,  "Unsubtracted"));
    p1.add_curve({smin, smax}, [&](double s){ return imag(pwave1->evaluate(s+IEPS) - pwave1->evaluate(0))/1E3; },                 dashed(jpacColor::Red));
    p1.add_curve({smin, smax}, [&](double s){ return real(pwave2->evaluate(s+IEPS) - pwave2->evaluate(0))/1E3; },                 solid( jpacColor::Green,"Once-subtracted"));
    p1.add_curve({smin, smax}, [&](double s){ return imag(pwave2->evaluate(s+IEPS) - pwave2->evaluate(0))/1E3; },                 dashed(jpacColor::Green));
    p1.shade_region({kin->sth(), kin->pth()});

    plot p2 = plotter.new_plot();
    p2.set_curve_points(400);
    p2.set_ranges({smin, smax}, {0, 5});
    p2.set_legend(0.68, 0.75);
    p2.set_labels("#sigma  [GeV^{2}]", "| #it{f}(#sigma) - #it{f}(0) | / 10^{3}");
    
    p2.add_curve({smin, smax}, [&](double s){ return abs(omnes_norm*pwave1->omnes(s+IEPS) - omnes_norm*pwave1->omnes(0))/1E3; }, solid(jpacColor::Blue, "Omnes"));
    p2.add_curve({smin, smax}, [&](double s){ return abs(pwave1->evaluate(s+IEPS) - pwave1->evaluate(0))/1E3; },                 solid(jpacColor::Red,  "Unsubtracted"));
    p2.add_curve({smin, smax}, [&](double s){ return abs(pwave2->evaluate(s+IEPS) - pwave2->evaluate(0))/1E3; },                 solid(jpacColor::Green,"Once-subtracted"));
    p2.shade_region({kin->sth(), kin->pth()});

    p1.save("reim_compare.pdf");
    p2.save("abs_compare.pdf");
    plotter.combine({2,1}, {p1,p2}, "compare.pdf");
   
};