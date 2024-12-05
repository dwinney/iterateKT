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
#include "decays/isohalf_pseudoscalar.hpp"
#include "plotter.hpp"

void kaon_decay()
{
    using namespace iterateKT;

    // -----------------------------------------------------------------------
    
    // Set up general kinematics so everything knows masses
    // Assume masses are given in terms of pion mass
    kinematics kin = new_kinematics(M_KAON, M_PION);

    // Set up our amplitude 
    amplitude amp = new_amplitude<charged_mode>(kin);

    // Add all the isobars, note the order they are added will be the order
    // the basis functions are generated
    amp->add_isobar<dI1_tI0_P1>(2); 

    // -----------------------------------------------------------------------
    // Iterate N times

    int N = 5;
    
    timer timer;
    timer.start();
    for (int i = 1; i <= N; i++)
    {
        amp->iterate();
        timer.lap("iteration " + std::to_string(i));
    }
    timer.stop();
    timer.print_elapsed();

    // -----------------------------------------------------------------------
    // Import data to 

    auto st1 = import_data<3>("scripts/bachir/St1.dat");

    // -----------------------------------------------------------------------
    // Plot Results

    plotter plotter;
    double smin =  +0.06;
    double smax =  +0.15;

    std::array<std::string,7> labels = {"0th", "1st", "2nd", "3rd", "4th", "5th", "6th"};

    auto plot_basis = [&](isobar isobar, int i, std::string label)
    {
        plot p = plotter.new_plot();
        p.set_curve_points(1000);
        p.set_legend(false);
        p.set_labels("#it{s} [GeV^{2}]", label);
       
        p.add_curve({smin, smax}, [&](double s){ return     std::real(isobar->basis_function(i, s+IEPS)); }, solid(jpacColor::Red));
        p.add_curve({smin, smax}, [&](double s){ return 100*std::imag(isobar->basis_function(i, s+IEPS)); }, solid(jpacColor::Blue));
        p.add_curve(st1[0], st1[1], dashed(jpacColor::DarkGrey));
        p.add_curve(st1[0], 100*st1[2], dashed(jpacColor::DarkGrey));
        return p;
    };

    // Grab out isobars for plotting
    isobar dI12_M1t = amp->get_isobar(id::dI1_tI0_P1);

    plot ft1 = plot_basis(dI12_M1t, 1, "#tilde{F}_{1}(s)");
    ft1.set_ranges({smin, smax}, {-0.02, 0.2});
    ft1.save("kaon_isobars.pdf");
};