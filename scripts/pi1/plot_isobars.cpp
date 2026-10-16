// Plot the solutions of the KT equations for different production mechinisms 
// which will contribute to the π₁ → 3π Daltiz plot
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Instituto de Ciencias Nucleares (ICN)
//               Universidad Nacional Autonoma de Mexico (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "plotter.hpp"
#include "settings.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"

void plot_isobars()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;
    
    // -----------------------------------------------------------------------
    // Operating options

    uint   N                =   10; // Number of iterations
    uint   asymptotic_power =    1; // sigma^N behavior of the isobar
    
    auto constant = [&](complex sigma){return 1.;};
    auto deck = [](double m3pi)
    {
        double t = -0.12; // GeV^2
        return [t,m3pi](complex sigma){return pi1::deck(t, m3pi*m3pi, sigma);};  
    };

    // -----------------------------------------------------------------------
    // Set up, shouldnt need to change anything below this line
    
    timer timer;
    timer.start();

    std::vector<amplitude> amps;

    amps.emplace_back(new_amplitude<pi1>(new_kinematics(1.3, M_PION)));
    amps.emplace_back(new_amplitude<pi1>(new_kinematics(1.6, M_PION)));
    amps.emplace_back(new_amplitude<pi1>(new_kinematics(1.9, M_PION)));
    amps.emplace_back(new_amplitude<pi1>(new_kinematics(2.2, M_PION)));

    for (auto amp : amps)
    {
        double m3pi = amp->get_kinematics()->M();
        std::string name = "#it{m}_{3#pi} = " + to_string(m3pi) + " GeV";
        isobar pwave = amp->add_isobar<P_wave>({ constant, deck(m3pi)}, asymptotic_power, id::P_wave, name);
        amp->iterate(N);
        timer.lap(name + " iterated!");
    };
    
    // -----------------------------------------------------------------------
    // Plot the 0th, 1st, and last iterations

    plotter plotter;

    std::array<double,2> bounds = {0., 1.5};

    plot p1 = plotter.new_plot();
    p1.set_curve_points(1000);
    p1.set_xrange(bounds);
    p1.set_labels("#sigma   [GeV^{2}]", "#it{F}#kern[-0.3]{_{1}} (#sigma + #it{i}#epsilon)");
    p1.set_legend(0.65, 0.65);

    isobar uniterated = amps[0]->get_isobar(id::P_wave);
    p1.add_curve(bounds, [&](double s) { return std::real(uniterated->omnes(s+IEPS)); }, solid(jpacColor::DarkGrey, "#Omega(#sigma)"));
    p1.add_curve(bounds, [&](double s) { return std::imag(uniterated->omnes(s+IEPS)); }, dashed(jpacColor::DarkGrey));

    for (auto amp : amps)
    {
        isobar pwave = amp->get_isobar(id::P_wave);
        p1.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, pwave->name());
        p1.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });
    };
    p1.save("fc_comparison.pdf");

    plot p2 = plotter.new_plot();
    p2.set_curve_points(1000);
    p2.set_xrange(bounds);
    p2.add_header("#it{t} = #minus 0.12 GeV^{2}");
    p2.set_labels("#sigma   [GeV^{2}]", "#it{F}#kern[-0.3]{_{#Delta}} (#it{t}, #it{m}_{3#pi}^{2} #; #sigma + #it{i}#epsilon)");
    p2.set_legend(0.65, 0.6);

    p2.add_curve(bounds, [&](double s) { return std::real(uniterated->omnes(s+IEPS)); }, solid(jpacColor::DarkGrey, "#Omega(#sigma)"));
    p2.add_curve(bounds, [&](double s) { return std::imag(uniterated->omnes(s+IEPS)); }, dashed(jpacColor::DarkGrey));

    for (auto amp : amps)
    {
        isobar pwave = amp->get_isobar(id::P_wave);
        p2.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(1, s+IEPS)); }, pwave->name());
        p2.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(1, s+IEPS)); });
    };
    p2.save("fd_comparison.pdf");

    timer.stop();
    timer.print_elapsed();
};