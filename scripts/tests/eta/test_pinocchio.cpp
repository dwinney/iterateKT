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
    amplitude amp_I1 = new_amplitude<I1_transition>(eta);
    amp_I1->add_isobar<I1_transition::S0_wave>(2);
    amp_I1->add_isobar<I1_transition::P1_wave>(1);
    amp_I1->add_isobar<I1_transition::S2_wave>(0);

    std::vector<isobar> previous = amp_I1->get_isobars();
    
    // Isolate our pwave
    isobar S0 = amp_I1->get_isobar(kI1_S0);
    isobar P1 = amp_I1->get_isobar(kI1_P1);
    isobar S2 = amp_I1->get_isobar(kI1_S2);

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    double smax = 35;

    plot p1 = plotter.new_plot();
    p1.set_curve_points(400);
    p1.set_legend(0.3, 0.4);
    p1.set_ranges({A, smax}, {-10, 6});
    p1.add_curve({A, smax}, [&](double s){ return std::imag(S0->pinocchio_integral(0, s, previous)); }, "Imaginary");
    p1.add_curve({A, smax}, [&](double s){ return std::real(S0->pinocchio_integral(0, s, previous)); }, "Real");
    p1.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{0}^{#alpha}(s)");

    plot p2 = plotter.new_plot();
    p2.set_curve_points(400);
    p2.set_ranges({A, smax}, {-20, 10});
    p2.add_curve({A, smax}, [&](double s){ return std::imag(P1->pinocchio_integral(0, s, previous)); });
    p2.add_curve({A, smax}, [&](double s){ return std::real(P1->pinocchio_integral(0, s, previous)); });
    p2.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{1}^{#alpha}(s)");

    plot p3 = plotter.new_plot();
    p3.set_curve_points(400);
    p3.set_ranges({A, smax}, {-15, 10});
    p3.add_curve({A, smax}, [&](double s){ return std::imag(S2->pinocchio_integral(0, s, previous)); });
    p3.add_curve({A, smax}, [&](double s){ return std::real(S2->pinocchio_integral(0, s, previous)); });
    p3.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{2}^{#alpha}(s)");

    plot p4 = plotter.new_plot();
    p4.set_curve_points(400);
    p4.set_legend(0.3, 0.4);
    p4.set_ranges({A, smax}, {-13, 40});
    p4.add_curve({A, smax}, [&](double s){ return std::imag(S0->pinocchio_integral(1, s, previous)); });
    p4.add_curve({A, smax}, [&](double s){ return std::real(S0->pinocchio_integral(1, s, previous)); });
    p4.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{0}^{#beta}(s)");

    plot p5 = plotter.new_plot();
    p5.set_curve_points(400);
    p5.set_ranges({A, smax}, {-200, 80});
    p5.add_curve({A, smax}, [&](double s){ return std::imag(P1->pinocchio_integral(1, s, previous)); });
    p5.add_curve({A, smax}, [&](double s){ return std::real(P1->pinocchio_integral(1, s, previous)); });
    p5.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{1}^{#beta}(s)");

    plot p6 = plotter.new_plot();
    p6.set_curve_points(400);
    p6.set_ranges({A, smax}, {-15, 50});
    p6.add_curve({A, smax}, [&](double s){ return std::imag(S2->pinocchio_integral(1, s, previous)); });
    p6.add_curve({A, smax}, [&](double s){ return std::real(S2->pinocchio_integral(1, s, previous)); });
    p6.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{2}^{#beta}(s)");

    plot p7 = plotter.new_plot();
    p7.set_curve_points(400);
    p7.set_legend(0.3, 0.4);
    p7.set_ranges({A, smax}, {-200, 190});
    p7.add_curve({A, smax}, [&](double s){ return std::imag(S0->pinocchio_integral(2, s, previous)); });
    p7.add_curve({A, smax}, [&](double s){ return std::real(S0->pinocchio_integral(2, s, previous)); });
    p7.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{0}^{#gamma}(s)");

    plot p8 = plotter.new_plot();
    p8.set_curve_points(400);
    p8.set_ranges({A, smax}, {-290, 90});
    p8.add_curve({A, smax}, [&](double s){ return std::imag(P1->pinocchio_integral(2, s, previous)); });
    p8.add_curve({A, smax}, [&](double s){ return std::real(P1->pinocchio_integral(2, s, previous)); });
    p8.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{1}^{#gamma}(s)");

    plot p9 = plotter.new_plot();
    p9.set_curve_points(400);
    p9.set_ranges({A, smax}, {-150, 200});
    p9.add_curve({A, smax}, [&](double s){ return std::imag(S2->pinocchio_integral(2, s, previous)); });
    p9.add_curve({A, smax}, [&](double s){ return std::real(S2->pinocchio_integral(2, s, previous)); });
    p9.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{2}^{#gamma}(s)");

    plotter.combine({3,3},{p1,p2,p3,p4,p5,p6,p7,p8,p9}, "pinocchio.pdf");

    timer.stop();
    timer.print_elapsed();
};