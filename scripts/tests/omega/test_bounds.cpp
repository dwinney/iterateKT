// Make sure the correct continuation of integration bounds is used
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
#include "decays/vector.hpp"

#include "plotter.hpp"

void test_bounds()
{
    using namespace iterateKT;

    // Set up general kinematics so everything knows masses
    
    kinematics omega = new_kinematics(M_OMEGA/M_PION, 1);
    double sth = omega->sth();
    double pth = omega->pth();
    double rth = omega->rth();
    double smax = 60;

    // -----------------------------------------------------------------------
    plotter plotter;

    plot p1 = plotter.new_plot();
    p1.set_curve_points(200);
    p1.set_legend(0.7, 0.6);
    p1.add_curve({sth, smax},[&](double s){ return std::real(omega->t_plus(s));},  "#it{t}_{#plus}(#it{s})");
    p1.add_curve({sth, smax},[&](double s){ return std::real(omega->t_minus(s));}, "#it{t}_{#minus}(#it{s})");
    p1.add_curve({pth, rth}, [&](double s){ return std::real(omega->t_curve(omega->phi_plus(s)));},  dashed(jpacColor::Orange, "#it{t}(#phi_{#plus}(#it{s}))"));
    p1.add_curve({pth, rth}, [&](double s){ return std::real(omega->t_curve(omega->phi_minus(s)));}, dotted(jpacColor::Green,  "#it{t}(#phi_{#minus}(#it{s}))"));
    p1.set_labels("#it{s} / m_{#pi}^{2}", "Re #it{t} / m_{#pi}^{2}");

    plot p2 = plotter.new_plot();
    p2.set_curve_points(200);
    p2.add_curve({sth, smax},[&](double s){ return std::imag(omega->t_plus(s));});
    p2.add_curve({sth, smax},[&](double s){ return std::imag(omega->t_minus(s));});
    p2.add_curve({pth, rth}, [&](double s){ return std::imag(omega->t_curve(omega->phi_plus(s)));},  dashed(jpacColor::Orange));
    p2.add_curve({pth, rth}, [&](double s){ return std::imag(omega->t_curve(omega->phi_minus(s)));}, dotted(jpacColor::Green));
    p2.set_labels("#it{s} / m_{#pi}^{2}", "Im #it{t} / m_{#pi}^{2}");

    plotter.combine({2,1}, {p1,p2}, "bounds.pdf");
};