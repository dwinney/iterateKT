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

    uint   N                =    9; // Number of iterations
    uint   asymptotic_power =    3; // sigma^N behavior of the isobar
    
    auto constant = [&](complex sigma){return 1.;};
    auto gen_deck = [](double m3pi)
    {
        double t = -0.1; // GeV^2
        return [t,m3pi](complex sigma){return pi1::deck(t, m3pi*m3pi, sigma);};  
    };
    auto gen_bubble = [](double m3pi, double lam)
    {
        return [m3pi,lam](complex sigma){return pi1::bubble(m3pi*m3pi, sigma, norm(lam));};  
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
        isobar pwave = amp->add_isobar<P_wave>({ constant, gen_bubble(m3pi, 0.05), gen_bubble(m3pi, 0.77), gen_deck(m3pi)}, asymptotic_power, id::P_wave, name);
        amp->iterate(N);
        timer.lap(name + " iterated!");
    };
    
    // -----------------------------------------------------------------------
    // Plot the 0th, 1st, and last iterations

    plotter plotter;


    std::array<double,2> bounds = {0., 2.0};

    plot p1 = plotter.new_plot();
    p1.set_curve_points(1000);
    p1.set_xrange(bounds);
    p1.set_labels("#sigma   [GeV^{2}]", "#it{F}#kern[-0.3]{_{1}} (#sigma + #it{i}#epsilon)");
    p1.set_legend(0.65, 0.65);

    isobar uniterated = amps[0]->get_isobar(id::P_wave);
    p1.add_curve(bounds, [&](double s) { return std::real(uniterated->basis_function(0, 0, s+IEPS)); }, solid(jpacColor::DarkGrey, "#Omega(#sigma)"));
    p1.add_curve(bounds, [&](double s) { return std::imag(uniterated->basis_function(0, 0, s+IEPS)); }, dashed(jpacColor::DarkGrey));

    for (auto amp : amps)
    {
        isobar pwave = amp->get_isobar(id::P_wave);
        p1.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, pwave->name());
        p1.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });
    };
    p1.save("f1_comparison.pdf");

    plot p2 = plotter.new_plot();
    p2.set_curve_points(1000);
    p2.set_xrange(bounds);
    p2.add_header("#it{t} = #minus 0.1 GeV^{2}");
    p2.set_labels("#sigma   [GeV^{2}]", "#it{F}#kern[-0.3]{_{#Delta}} (#it{t}, #it{m}_{3#pi}^{2} #; #sigma + #it{i}#epsilon)");
    p2.set_legend(0.6, 0.6);

    p2.add_curve(bounds, [&](double s) { return std::real(uniterated->basis_function(0, 3, s+IEPS)); }, solid(jpacColor::DarkGrey, "#Delta(t, #it{m}_{3#pi}; #sigma) #Omega(#sigma)"));
    p2.add_curve(bounds, [&](double s) { return std::imag(uniterated->basis_function(0, 3, s+IEPS)); }, dashed(jpacColor::DarkGrey));

    for (auto amp : amps)
    {
        isobar pwave = amp->get_isobar(id::P_wave);
        p2.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(3, s+IEPS)); }, pwave->name());
        p2.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(3, s+IEPS)); });
    };
    p2.save("fd_comparison.pdf");

    plot p3 = plotter.new_plot();
    p3.set_curve_points(1000);
    p3.set_xrange(bounds);
    p3.add_header("#it{m}_{3#pi} = 1.6 GeV");
    p3.set_labels("#sigma   [GeV^{2}]", "#it{F}#kern[-0.3]{_{#it{B}}} ({#it{s}, #it{t}}#; #sigma + #it{i}#epsilon) / #it{F}#kern[-0.3]{_{#it{B}}} ({#it{s}, #it{t}}#; 0)");
    p3.set_legend(0.55, 0.6);

    std::array<std::string,4> labels = {"Contact", "Bubble (#Lambda = 200 MeV)", "Bubble (#Lambda = 770 MeV)", "Deck (#it{t} = #minus 0.1 GeV^{2})"};
    for (int i = 0; i < 4; i++)
    {
        isobar pwave = amps[1]->get_isobar(id::P_wave);
        complex norm = pwave->basis_function(i, 0);
        p3.add_curve( bounds, [&](double s) { return std::real(pwave->basis_function(i, s+IEPS)/norm); }, labels[i]);
        p3.add_dashed(bounds, [&](double s) { return std::imag(pwave->basis_function(i, s+IEPS)/norm); });
    };
    p3.save("fb_comparison.pdf");


    timer.stop();
    timer.print_elapsed();
};