// Evaluate the pinocchio integration
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


void test_F1_inhomogeneity()
{
    using namespace iterateKT;
    using namespace pseudoscalar;

    // Set up general kinematics so everything knows masses
    // Assume masses are given in terms of pion mass
    kinematics kin = new_kinematics(M_ETA/M_PION, 1.);

    // Significant points in integration path
    double A = kin->A();
    double B = kin->B();
    double C = kin->C();
    double D = kin->D();

    // Set up our amplitude 
    amplitude amp_I1 = new_amplitude<I1_transition>(kin);
    amp_I1->add_isobar<I1_transition::S0_wave>(2);
    amp_I1->add_isobar<I1_transition::P1_wave>(1);
    amp_I1->add_isobar<I1_transition::S2_wave>(0);

    std::vector<isobar> previous = amp_I1->get_isobars();
    
    // Isolate our pwave
    isobar P1 = amp_I1->get_isobar(kI1_P1);

    // -----------------------------------------------------------------------
    
    timer timer;
    plotter plotter;

    timer.start();

    auto grid = P1->calculate_next(previous);
    timer.lap("grid");

    // first nontrivial iteration
    iteration first = new_iteration(kin, grid, pseudoscalar::default_settings());

    double smax = 20;

    plot p1 = plotter.new_plot();
    p1.set_legend(0.25, 0.2);
    p1.set_curve_points(100);
    p1.add_curve( {A, smax}, [&](double s){ return std::real(first->half_regularized_integrand(0, s)); }, solid(jpacColor::Blue, "Real"));
    p1.add_curve( {A, smax}, [&](double s){ return std::imag(first->half_regularized_integrand(0, s)); }, solid(jpacColor::Red,  "Imaginary"));
    p1.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{1}^{#alpha,0} / #nu^{3}");

    plot p2 = plotter.new_plot();
    p2.set_legend(0.25, 0.2);
    p2.set_curve_points(300);
    p2.add_curve( {A, smax}, [&](double s){ return std::real(first->regularized_integrand(0, s)); }, solid(jpacColor::Green,  "Real"));
    p2.add_curve( {A, smax}, [&](double s){ return std::imag(first->regularized_integrand(0, s)); }, solid(jpacColor::Orange, "Imaginary"));
    p2.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{1}^{#alpha,0} / #kappa^{3}");

    // plotter.combine({2,1}, {p1,p2}, "F1_inhomogeneity.pdf");

    plot p3 = plotter.new_plot();
    p3.set_legend(0.25, 0.2);
    p3.set_curve_points(300);
    p3.add_curve( {A, smax}, [&](double s){ return std::real(first->integral(0, s+IEPS)); }, solid(jpacColor::Green,  "Real"));
    p3.add_curve( {A, smax}, [&](double s){ return std::imag(first->integral(0, s+IEPS)); }, solid(jpacColor::Orange, "Imaginary"));
    p3.set_labels("#it{s} / #it{m}_{#pi}^{2}", "I_{1}^{#alpha,1}");

    plotter.combine({3,1}, {p1,p2,p3}, "F1_inhomogeneity.pdf");

    // plot p3 = plotter.new_plot();
    // p3.set_legend(0.45, 0.6);
    // p3.set_curve_points(100);
    // p3.add_curve( {A, smax}, [&](double s){ return std::real(first->half_regularized_integrand(1, s)); }, solid(jpacColor::Blue));
    // p3.add_curve( {A, smax}, [&](double s){ return std::imag(first->half_regularized_integrand(1, s)); }, solid(jpacColor::Red));
    // p3.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{1}^{#beta,0} / #nu^{3}");

    // plot p4 = plotter.new_plot();
    // p4.set_legend(0.45, 0.4);
    // p4.set_curve_points(100);
    // p4.add_curve( {A, smax}, [&](double s){ return std::real(first->regularized_integrand(1, s)); }, solid(jpacColor::Green));
    // p4.add_curve( {A, smax}, [&](double s){ return std::imag(first->regularized_integrand(1, s)); }, solid(jpacColor::Orange));
    // p4.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{1}^{#beta,0} / #kappa^{3}");

    // plot p5 = plotter.new_plot();
    // p5.set_legend(0.45, 0.6);
    // p5.set_curve_points(100);
    // p5.add_curve( {A, smax}, [&](double s){ return std::real(first->half_regularized_integrand(1, s)); }, solid(jpacColor::Blue));
    // p5.add_curve( {A, smax}, [&](double s){ return std::imag(first->half_regularized_integrand(1, s)); }, solid(jpacColor::Red));
    // p5.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{1}^{#gamma,0} / #nu^{3}");

    // plot p6 = plotter.new_plot();
    // p6.set_legend(0.45, 0.4);
    // p6.set_curve_points(100);
    // p6.add_curve( {A, smax}, [&](double s){ return std::real(first->regularized_integrand(1, s)); }, solid(jpacColor::Green));
    // p6.add_curve( {A, smax}, [&](double s){ return std::imag(first->regularized_integrand(1, s)); }, solid(jpacColor::Orange));
    // p6.set_labels("#it{s} / #it{m}_{#pi}^{2}", "#tilde{F}_{1}^{#gamma,0} / #kappa^{3}");

    // plotter.combine({2,3}, {p1,p2,p3,p4,p5,p6}, "F1_inhomogeneity.pdf");

    timer.stop();
    timer.print_elapsed();
};