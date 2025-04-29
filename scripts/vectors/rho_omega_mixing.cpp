// Twice subtracted KT amplitudes for omega decay with only P-wave in [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2006.01058
// ------------------------------------------------------------------------------


#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"

#include "isobars/vector.hpp"
#include "amplitudes/vector.hpp"

#include "plotter.hpp"

void rho_omega_mixing()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // Set up general kinematics so everything knows masses
    // Use masses in units of GeV
    double M_V = 3.096;
    kinematics kinematics = new_kinematics(M_V, M_PION);
    
    // Significant points in integration path
    double sth = kinematics->sth(); 
    double pth = kinematics->pth(); 

    // Set up our amplitude 
    amplitude amplitude = new_amplitude<iterateKT::rho_omega_mixing>(kinematics);

    // Settings to adjust threshold from the defaults which are set for the omega
    settings jpsi_settings = P_wave::default_settings();
    jpsi_settings._intermediate_energy  = 15.;
    jpsi_settings._cutoff               = 40.;
    jpsi_settings._interpolation_points = {400, 10, 300};
    jpsi_settings._angular_integrator_depth = 6;

    auto constant = [](complex s){ return complex(1.); };

    // Charged channel is just a normal unsubtracted dispersion relation
    amplitude->add_isobar<charged>(constant, 1, id::charged, "charged", jpsi_settings);

    // The neutral one however takes the omega BW in addition to the constant
    auto omega_bw = [](complex s)
    {
        double m = M_OMEGA, g = 8.68E-3;
        return s/(s-m*m+I*m*g);
    };
    jpsi_settings._angular_integrator_depth = 0;
    amplitude->add_isobar<neutral>({constant, omega_bw}, 1, id::neutral, "neutral", jpsi_settings);

    isobar charged = amplitude->get_isobar(id::charged);
    isobar neutral = amplitude->get_isobar(id::neutral);

    // -----------------------------------------------------------------------
    bool do_next           = true;
    std::string next_label = "1st";
    std::vector<std::array<std::string,2>> previous_iterations = {
        // {"1st_charged.dat", "1st_neutral.dat"},
        // {"2nd_charged.dat", "2nd_neutral.dat"},
        // {"3rd_charged.dat", "3rd_neutral.dat"},
        // {"4th_charged.dat", "4th_neutral.dat"},
    };
    timer timer;
    timer.start();
    for (auto file : previous_iterations)
    {
        std::string path = "scripts/vectors/";
        charged->import_iteration<3>(path+file[0]);
        neutral->import_iteration<3>(path+file[1]);
    };
    timer.lap("Imported");
    if (do_next)
    {
        amplitude->iterate();
        timer.lap("Iterated");
        amplitude->export_solution(next_label);
        timer.lap("Exported");
    }

    // -----------------------------------------------------------------------

    plotter plotter;
    
    std::array<double,2> bounds = {0., 1};
    plot p1 = plotter.new_plot();
    p1.set_legend(0.25, 0.7);
    p1.set_curve_points(1000);
    p1.set_labels("#sigma [GeV^{2}]", "f_{#alpha}(#sigma + #it{i}#epsilon)");
    
    plot p2 = plotter.new_plot();
    p2.set_legend(0.25, 0.7);
    p2.set_curve_points(1000);
    p2.set_labels("#sigma [GeV^{2}]", "f_{#beta}(#sigma + #it{i}#epsilon)");
    
    plot p3 = plotter.new_plot();
    p3.set_legend(0.25, 0.7);
    p3.set_curve_points(1000);
    p3.set_labels("#sigma [GeV^{2}]", "f_{#gamma}(#sigma + #it{i}#epsilon)");
    
    plot p4 = plotter.new_plot();
    p4.set_legend(0.25, 0.7);
    p4.set_curve_points(2000);
    p4.set_labels("#sigma [GeV^{2}]", "g_{#alpha}(#sigma + #it{i}#epsilon)");
    
    plot p5 = plotter.new_plot();
    p5.set_legend(0.25, 0.7);
    p5.set_curve_points(2000);
    p5.set_labels("#sigma [GeV^{2}]", "g_{#beta}(#sigma + #it{i}#epsilon)");
    
    plot p6 = plotter.new_plot();
    p6.set_legend(0.25, 0.7);
    p6.set_curve_points(2000);
    p6.set_labels("#sigma [GeV^{2}]", "g_{#gamma}(#sigma + #it{i}#epsilon)");

    std::array<std::string,5> labels = {"0th", "1st", "2nd", "3rd", "4th"};
    for (int i = 0; i <= previous_iterations.size()+do_next; i++)
    {
        p1.add_curve( bounds, [&](double s) { return std::real(charged->basis_function(i, 0, s+IEPS)); }, labels[i]);
        p1.add_dashed(bounds, [&](double s) { return std::imag(charged->basis_function(i, 0, s+IEPS)); });
        p2.add_curve( bounds, [&](double s) { return std::real(charged->basis_function(i, 1, s+IEPS)); }, labels[i]);
        p2.add_dashed(bounds, [&](double s) { return std::imag(charged->basis_function(i, 1, s+IEPS)); });
        p3.add_curve( bounds, [&](double s) { return std::real(charged->basis_function(i, 2, s+IEPS)); }, labels[i]);
        p3.add_dashed(bounds, [&](double s) { return std::imag(charged->basis_function(i, 2, s+IEPS)); });
        
        p4.add_curve( bounds, [&](double s) { return std::real(neutral->basis_function(i, 0, s+IEPS)); }, labels[i]);
        p4.add_dashed(bounds, [&](double s) { return std::imag(neutral->basis_function(i, 0, s+IEPS)); });
        p5.add_curve( bounds, [&](double s) { return std::real(neutral->basis_function(i, 1, s+IEPS)); }, labels[i]);
        p5.add_dashed(bounds, [&](double s) { return std::imag(neutral->basis_function(i, 1, s+IEPS)); });
        p6.add_curve( bounds, [&](double s) { return std::real(neutral->basis_function(i, 2, s+IEPS)); }, labels[i]);
        p6.add_dashed(bounds, [&](double s) { return std::imag(neutral->basis_function(i, 2, s+IEPS)); });
    };
    
    plotter.combine({3,2}, {p1,p2,p3,p4,p5,p6}, "rho_omega_isobars.pdf");  
    timer.lap("Plotted");

    timer.stop();
    timer.print_elapsed();

};