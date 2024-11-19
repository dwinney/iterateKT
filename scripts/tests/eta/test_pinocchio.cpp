// Evaluate the pinocchio integration for C-conserving eta -> 3pi
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "basis.hpp"
#include "decays/pseudoscalar.hpp"

#include "plotter.hpp"


void test_pinocchio()
{
    using namespace iterateKT;
    using namespace pseudoscalar;

    // Set up general kinematics so everything knows masses
    // Assume masses are given in terms of pion mass
    kinematics eta = new_kinematics(M_ETA/M_PION, 1.);

    // Significant points in integration path
    double A = eta->A();
    double B = eta->B();
    double C = eta->C();
    double D = eta->D();

    // Set up our amplitude 
    amplitude amp_dI1 = new_amplitude<neutral_mode>(eta);
    amp_dI1->add_isobar<dI1_S0>(2);
    amp_dI1->add_isobar<dI1_P1>(1);
    amp_dI1->add_isobar<dI1_S2>(0);


    
    // Isolate our pwave
    isobar S0 = amp_dI1->get_isobar(kdI1_S0);
    isobar P1 = amp_dI1->get_isobar(kdI1_P1);
    isobar S2 = amp_dI1->get_isobar(kdI1_S2);

    // -----------------------------------------------------------------------
    
    plotter plotter;

    timer timer;

    timer.start();
    amp_dI1->iterate();
    timer.lap("Iterated");

    std::vector<isobar> previous = amp_dI1->get_isobars();

    double smax = 35;

    auto plot_pinocchio = [&](isobar isobar, int i, std::string label)
    {
        plot p = plotter.new_plot();    
        p.set_curve_points(100);
        p.print_to_terminal(true);
        p.set_legend(false);
        p.add_curve( {A, smax}, [&](double s){ return std::imag(isobar->pinocchio_integral(i, s, previous)); }, "Imaginary");
        p.add_curve( {A, smax}, [&](double s){ return std::real(isobar->pinocchio_integral(i, s, previous)); }, "Real");
        p.set_labels("#it{s} / #it{m}_{#pi}^{2}", label);
        timer.lap("plot "+std::to_string(i));
        return p;
    };

    plot f0a = plot_pinocchio(S0, 0, "#tilde{F}_{0}^{#alpha}(s)");
    plot f1a = plot_pinocchio(P1, 0, "#tilde{F}_{1}^{#alpha}(s)");
    plot f2a = plot_pinocchio(S2, 0, "#tilde{F}_{0}^{#alpha}(s)");
    plot f0b = plot_pinocchio(S0, 1, "#tilde{F}_{0}^{#beta}(s)");
    plot f1b = plot_pinocchio(P1, 1, "#tilde{F}_{1}^{#beta}(s)");
    plot f2b = plot_pinocchio(S2, 1, "#tilde{F}_{2}^{#beta}(s)");
    plot f0g = plot_pinocchio(S0, 2, "#tilde{F}_{0}^{#gamma}(s)");
    plot f1g = plot_pinocchio(P1, 2, "#tilde{F}_{1}^{#gamma}(s)");
    plot f2g = plot_pinocchio(S2, 2, "#tilde{F}_{2}^{#gamma}(s)");

    plotter.combine({3,3}, {f0a, f1a, f2a, f0b, f1b, f2b, f0g, f1g, f2g}, "pinocchio.pdf");

    timer.stop();
    timer.print_elapsed();
};