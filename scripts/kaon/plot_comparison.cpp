// Basis functions for K -> 3pi following [1]
// These are calcualted in the isospin limit.
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

#include <filesystem>
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "basis.hpp"
#include "plotter.hpp"

#include "amplitudes/kaon.hpp"
#include "isobars/kaon.hpp"

void plot_comparison()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // -----------------------------------------------------------------------
    
    // Isospin averaged Kaon and charged pion
    kinematics kin = new_kinematics(M_KAON_AVG, M_PION_PM);

    // Set up our amplitude 
    solver solver(kin);

    // Add all the isobars, note the order they are added will be the order
    // the basis functions are generated
    std::vector<uint> empty = {};
    isobar F0 = solver.add_isobar<I1_S0>({0, 1, 2}, 2, id::dI1_I1_S0, "F0");
    isobar F1 = solver.add_isobar<I1_P1>({1},       1, id::dI1_I1_P1, "F1");
    isobar F2 = solver.add_isobar<I1_S2>(empty,     2, id::dI1_I1_S2, "F2"); 
    isobar G1 = solver.add_isobar<I0_P1>({1},       1, id::dI1_I0_P1, "G1");
    isobar H1 = solver.add_isobar<I2_P1>({0, 1},    1, id::dI3_I2_P1, "H1");
    isobar H2 = solver.add_isobar<I2_S2>(empty,     2, id::dI3_I2_S2, "H2");

    // -----------------------------------------------------------------------
    // Also grab the isobars from [1]

    auto S0  = import_data<7>("data/K_3pi/orsay_isobars/S_0.dat");
    auto S1  = import_data<7>("data/K_3pi/orsay_isobars/S_1.dat");
    auto S2  = import_data<7>("data/K_3pi/orsay_isobars/S_2.dat");
    auto S3  = import_data<7>("data/K_3pi/orsay_isobars/S_3.dat");
    auto Sp  = import_data<3>("data/K_3pi/orsay_isobars/St1.dat");
    auto St1 = import_data<5>("data/K_3pi/orsay_isobars/StI_0.dat");
    auto St2 = import_data<5>("data/K_3pi/orsay_isobars/StI_1.dat");

    // -----------------------------------------------------------------------
    // Plot Results

    // Path to precalculated isobar files
    std::string path   = "/scripts/kaon/basis_functions/";
    std::string prefix = "basis_";
    // Import everything 
    for (auto iso : solver.get_isobars()) iso->import_iteration<7>(path+prefix+iso->name()+".dat");

    plotter plotter;
    double smin =  +0.05;
    double smax =  +0.15;

    auto plot_basis = [&](isobar isobar, int i, std::string func, std::array<double,2> ranges = {0, 0})
    {
        plot p = plotter.new_plot();
        p.set_curve_points(1000);
        p.set_legend(false);
        (ranges[0] == ranges[1]) ? p.set_xrange({smin, smax}) : p.set_ranges({smin, smax}, ranges);

        auto rF = [&](){return [isobar,i](double s){return real(isobar->basis_function(i, s+IEPS));}; };
        auto iF = [&](){return [isobar,i](double s){return imag(isobar->basis_function(i, s+IEPS));}; };

        std::string label = func+"(#it{s} + #it{i}#epsilon)";
        p.set_labels("#it{s}  [GeV^{2}]", label);
        p.shade_region({kin->sth(), kin->pth()});
        p.add_vertical(kin->rth());
        p.add_curve({smin, smax}, rF(), solid(jpacColor::Red));
        p.add_curve({smin, smax}, iF(), solid(jpacColor::Blue));
        return p;
    };

    std::vector<plot> F_plots;
    F_plots.emplace_back(plot_basis(F0, 0, "F_{0}^{#alpha}", {-0.1, 1.7}));
    F_plots.back().add_curve(S0[0], S0[1], dashed(jpacColor::DarkGrey));
    F_plots.back().add_curve(S0[0], S0[2], dashed(jpacColor::DarkGrey));
    F_plots.emplace_back(plot_basis(F1, 0, "F_{1}^{#alpha}", {-50E-3, 150E-3}));
    F_plots.back().add_curve(S0[0], S0[3], dashed(jpacColor::DarkGrey));
    F_plots.back().add_curve(S0[0], S0[4], dashed(jpacColor::DarkGrey));
    F_plots.emplace_back(plot_basis(F2, 0, "F_{2}^{#alpha}", {-130E-3, 25E-3}));
    F_plots.back().add_curve(S0[0], S0[5], dashed(jpacColor::DarkGrey));
    F_plots.back().add_curve(S0[0], S0[6], dashed(jpacColor::DarkGrey));

    F_plots.emplace_back(plot_basis(F0, 1, "F_{0}^{#beta}",  {-20E-3, 275E-3}));
    F_plots.back().add_curve(S1[0], S1[1], dashed(jpacColor::DarkGrey));
    F_plots.back().add_curve(S1[0], S1[2], dashed(jpacColor::DarkGrey));
    F_plots.emplace_back(plot_basis(F1, 1, "F_{1}^{#beta}",  {-20E-4, 300E-4}));
    F_plots.back().add_curve(S1[0], S1[3], dashed(jpacColor::DarkGrey));
    F_plots.back().add_curve(S1[0], S1[4], dashed(jpacColor::DarkGrey));
    F_plots.emplace_back(plot_basis(F2, 1, "F_{2}^{#beta}",  {-10E-3, 6E-3}));
    F_plots.back().add_curve(S1[0], S1[5], dashed(jpacColor::DarkGrey));
    F_plots.back().add_curve(S1[0], S1[6], dashed(jpacColor::DarkGrey));
    
    F_plots.emplace_back(plot_basis(F0, 2, "F_{0}^{#gamma}", {-3E-3, 40E-3}));
    F_plots.back().add_curve(S2[0], S2[1], dashed(jpacColor::DarkGrey));
    F_plots.back().add_curve(S2[0], S2[2], dashed(jpacColor::DarkGrey));
    F_plots.emplace_back(plot_basis(F1, 2, "F_{1}^{#gamma}", {-50E-4, 40E-4}));
    F_plots.back().add_curve(S2[0], S2[3], dashed(jpacColor::DarkGrey));
    F_plots.back().add_curve(S2[0], S2[4], dashed(jpacColor::DarkGrey));
    F_plots.emplace_back(plot_basis(F2, 2, "F_{2}^{#gamma}", {-65E-5, 25E-5}));
    F_plots.back().add_curve(S2[0], S2[5], dashed(jpacColor::DarkGrey));
    F_plots.back().add_curve(S2[0], S2[6], dashed(jpacColor::DarkGrey));

    F_plots.emplace_back(plot_basis(F0, 3, "F_{0}^{#zeta}", {-5E-4, 40E-4}));
    F_plots.back().add_curve(S3[0], S3[1], dashed(jpacColor::DarkGrey));
    F_plots.back().add_curve(S3[0], S3[2], dashed(jpacColor::DarkGrey));
    F_plots.emplace_back(plot_basis(F1, 3, "F_{1}^{#zeta}", {-10E-3, 225E-3}));
    F_plots.back().add_curve(S3[0], S3[3], dashed(jpacColor::DarkGrey));
    F_plots.back().add_curve(S3[0], S3[4], dashed(jpacColor::DarkGrey));
    F_plots.emplace_back(plot_basis(F2, 3, "F_{2}^{#zeta}", {-125E-5, 75E-5}));
    F_plots.back().add_curve(S3[0], S3[5], dashed(jpacColor::DarkGrey));
    F_plots.back().add_curve(S3[0], S3[6], dashed(jpacColor::DarkGrey));
    plotter.combine({3,4}, F_plots, "Fs.pdf");

    plot pGr = plotter.new_plot();
    pGr.set_curve_points(1000);
    pGr.set_labels("#it{s}  [GeV^{2}]", "Re G_{1}^{#eta}(#it{s} + #it{i}#epsilon)");
    pGr.shade_region({kin->sth(), kin->pth()});
    pGr.add_vertical(kin->rth());
    pGr.set_ranges({smin, smax}, {-1, 1});
    pGr.add_curve({smin, smax}, [&](double s){ return real(G1->basis_function(4, s+IEPS)); }, solid(jpacColor::Red));
    pGr.add_curve(Sp[0], Sp[1], dashed(jpacColor::DarkGrey));

    plot pGi = plotter.new_plot();
    pGi.set_curve_points(1000);
    pGi.set_labels("#it{s}  [GeV^{2}]", "Im G_{1}^{#eta}(#it{s} + #it{i}#epsilon)");
    pGi.shade_region({kin->sth(), kin->pth()});
    pGi.add_vertical(kin->rth());
    pGi.set_ranges({smin, smax}, {-1E-4, 3E-4});
    pGi.add_curve({smin, smax}, [&](double s){ return imag(G1->basis_function(4, s+IEPS)); }, solid(jpacColor::Blue));
    pGi.add_curve(Sp[0], Sp[2], dashed(jpacColor::DarkGrey));
    plotter.combine({2,1}, {pGr, pGi}, "G.pdf");

    std::vector<plot> H_plots;
    H_plots.emplace_back(plot_basis(H1, 5, "H_{1}^{#mu}", {-1E-3, 1.5}));
    H_plots.back().add_curve(St1[0], St1[1], dashed(jpacColor::DarkGrey));
    H_plots.back().add_curve(St1[0], St1[2], dashed(jpacColor::DarkGrey));
    H_plots.emplace_back(plot_basis(H2, 5, "H_{2}^{#mu}", {-5E-3, 4E-2}));
    H_plots.back().add_curve(St1[0], St1[3], dashed(jpacColor::DarkGrey));
    H_plots.back().add_curve(St1[0], St1[4], dashed(jpacColor::DarkGrey));
    H_plots.emplace_back(plot_basis(H1, 6, "H_{1}^{#nu}", {-5E-3, 25E-2}));
    H_plots.back().add_curve(St2[0], St2[1], dashed(jpacColor::DarkGrey));
    H_plots.back().add_curve(St2[0], St2[2], dashed(jpacColor::DarkGrey));
    H_plots.emplace_back(plot_basis(H2, 6, "H_{2}^{#nu}", {-5E-3, 5E-3}));
    H_plots.back().add_curve(St2[0], St2[3], dashed(jpacColor::DarkGrey));
    H_plots.back().add_curve(St2[0], St2[4], dashed(jpacColor::DarkGrey));
    plotter.combine({2,2}, H_plots, "Hs.pdf");
};