// Basis functions for K -> 3pi following [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
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
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "basis.hpp"
#include "plotter.hpp"
#include "solver.hpp"

#include "isobars/kaon.hpp"

void kaon_decay()
{
    using namespace iterateKT;

    // -----------------------------------------------------------------------
    
    // Set up general kinematics so everything knows masses
    // Assume masses are given in terms of pion mass
    kinematics kin = new_kinematics(M_KAON, M_PION);

    // Set up our amplitude 
    solver solver(kin);

    // Add all the isobars, note the order they are added will be the order
    // the basis functions are generated
    solver.add_isobar<dI1_tI1_S0>({0, 1, 2}, 2, id::dI1_tI1_S0); // M0
    solver.add_isobar<dI1_tI1_P1>({1},       1, id::dI1_tI1_P1); // M1
    solver.add_isobar<dI1_tI1_S2>({},        2, id::dI1_tI1_S2); // M2
    solver.add_isobar<dI3_tI2_P1>({0, 1},    1, id::dI3_tI2_P1); // Nt1
    solver.add_isobar<dI3_tI2_S2>({},        2, id::dI3_tI2_S2); // Nt2
    solver.add_isobar<dI1_tI0_P1>({1},       1, id::dI1_tI0_P1); // Mt1

    // -----------------------------------------------------------------------
    // Iterate N times

    int N = 4;
    
    timer timer;
    timer.start();
    for (int i = 1; i <= N; i++)
    {
        solver.iterate();
        timer.lap("iteration " + std::to_string(i));
    }
    timer.stop();
    timer.print_elapsed();

    // -----------------------------------------------------------------------
    // Plot Results

    plotter plotter;
    double smin =  +0.06;
    double smax =  +0.15;

    auto plot_basis = [&](isobar isobar, std::array<double,4> scale, std::array<std::string,4> labels)
    {
        plot p = plotter.new_plot();
        p.set_curve_points(1000);
        p.set_labels("#it{s} [GeV^{2}]", "");

        p.add_curve({smin, smax}, [&](double s){ return scale[0]*std::real(isobar->basis_function(0, s+IEPS)); }, solid(jpacColor::Brown, labels[0]));
        p.add_curve({smin, smax}, [&](double s){ return scale[0]*std::imag(isobar->basis_function(0, s+IEPS)); }, dashed(jpacColor::Brown));

        p.add_curve({smin, smax}, [&](double s){ return scale[1]*std::real(isobar->basis_function(1, s+IEPS)); }, solid(jpacColor::Blue, labels[1]));
        p.add_curve({smin, smax}, [&](double s){ return scale[1]*std::imag(isobar->basis_function(1, s+IEPS)); }, dashed(jpacColor::Blue));

        p.add_curve({smin, smax}, [&](double s){ return scale[2]*std::real(isobar->basis_function(2, s+IEPS)); }, solid(jpacColor::Green, labels[2]));
        p.add_curve({smin, smax}, [&](double s){ return scale[2]*std::imag(isobar->basis_function(2, s+IEPS)); }, dashed(jpacColor::Green));

        p.add_curve({smin, smax}, [&](double s){ return scale[3]*std::real(isobar->basis_function(3, s+IEPS)); }, solid(jpacColor::Pink, labels[3]));
        p.add_curve({smin, smax}, [&](double s){ return scale[3]*std::imag(isobar->basis_function(3, s+IEPS)); }, dashed(jpacColor::Pink));
        return p;
    };

    auto plot_basis_tilde = [&](isobar isobar, std::array<double,2> scale, std::array<std::string,2> labels)
    {
        plot p = plotter.new_plot();
        p.set_curve_points(1000);
        p.set_labels("#it{s} [GeV^{2}]", "");

        p.add_curve({smin, smax}, [&](double s){ return scale[0]*std::real(isobar->basis_function(4, s+IEPS)); }, solid(jpacColor::Brown, labels[0]));
        p.add_curve({smin, smax}, [&](double s){ return scale[0]*std::imag(isobar->basis_function(4, s+IEPS)); }, dashed(jpacColor::Brown));

        p.add_curve({smin, smax}, [&](double s){ return scale[1]*std::real(isobar->basis_function(5, s+IEPS)); }, solid(jpacColor::Blue, labels[1]));
        p.add_curve({smin, smax}, [&](double s){ return scale[1]*std::imag(isobar->basis_function(5, s+IEPS)); }, dashed(jpacColor::Blue));
        return p;
    };

    // Grab out isobars for plotting
    isobar M0  = solver.get_isobar(id::dI1_tI1_S0);
    isobar M1  = solver.get_isobar(id::dI1_tI1_P1);
    isobar M2  = solver.get_isobar(id::dI1_tI1_S2);
    isobar Nt1 = solver.get_isobar(id::dI3_tI2_P1);
    isobar Nt2 = solver.get_isobar(id::dI3_tI2_S2);
    isobar Mt1 = solver.get_isobar(id::dI1_tI0_P1);

    // Labels for the legends of all the plots
    std::array<std::string,4> f0_labels =  {"S_{0}^{0}", "10 S_{0}^{1}", "10^{2} S_{0}^{2}", "10^{3} S_{0}^{3}"};
    std::array<std::string,4> f1_labels =  {"S_{1}^{0}", "10 S_{1}^{1}", "10^{2} S_{1}^{2}", "S_{1}^{3}"};
    std::array<std::string,4> f2_labels =  {"S_{2}^{0}", "10 S_{2}^{1}", "10^{2} S_{2}^{2}", "S_{2}^{3}"} ;


    // S functions
    plot f0 = plot_basis(M0, {1, 10, 1E2, 1E3}, f0_labels);
    f0.set_ranges({smin, smax}, {-0.3, 3.0});
    f0.set_legend(0.25, 0.7);

    plot f1 = plot_basis(M1, {1, 10, 1E2, 1},  f1_labels);
    f1.set_ranges({smin, smax}, {-0.4, 0.3});
    f1.set_legend(0.25, 0.2);

    plot f2 = plot_basis(M2, {1, 10, 1E2, 1E2}, f2_labels);
    f2.set_ranges({smin, smax}, {-0.12, 0.06});
    f2.set_legend(0.25, 0.2);

    // Stilde functions
    plot f2t = plot_basis_tilde(Nt2, {1, 10}, {"#tilde{S}_{2}^{0}", "10 #tilde{S}_{2}^{1}"});
    f2t.set_legend(0.25, 0.7);
    f2t.set_ranges({smin,smax}, {-0.04, 0.05});
    
    plot f1t = plot_basis_tilde(Nt1, {1, 10}, {"#tilde{S}_{1}^{0}", "10 #tilde{S}_{1}^{1}"});
    f1t.set_legend(0.25, 0.7);
    f1t.set_ranges({smin,smax}, {-0.06, 2});

    plot s1t = plotter.new_plot();
    s1t.set_curve_points(1000);
    s1t.set_labels("#it{s} [GeV^{2}]", "");
    s1t.set_legend(0.25, 0.7);
    s1t.set_ranges({smin,smax}, {-0.01, 0.2});
    s1t.add_curve({smin, smax}, [&](double s){ return     std::real(Mt1->basis_function(6, s+IEPS)); }, solid(jpacColor::Blue,  "Re #tilde{S}_{1}"));
    s1t.add_curve({smin, smax}, [&](double s){ return 100*std::imag(Mt1->basis_function(6, s+IEPS)); }, dashed(jpacColor::Blue, "10^{2} Im #tilde{S}_{1}"));

    plotter.combine({3,2}, {f0,f1,f2,f2t,f1t,s1t}, "kaon_isobars.pdf");
};