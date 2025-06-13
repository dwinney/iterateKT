//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2212.11767
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "solver.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "plotter.hpp"
#include "settings.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"

void kt_deck()
{
    using namespace iterateKT;
    using iterateKT::complex;

    uint   N    = 4;    // Number of iterations
    double t    = -0.1; // Production t
    double m3pi = 1.40; // 3pi invariant mass

    kinematics kinematics = new_kinematics(m3pi, M_PION);
    solver solver(kinematics);

    // The projection function is given by our Deck loop 
    double Delta0 = imag(pi1::projected_deck(t, m3pi*m3pi, 0)); // normalize to sigma = 0
    auto   Delta  = [&](complex sigma){return imag(pi1::projected_deck(t, m3pi*m3pi, sigma))/Delta0;};

    settings sets = default_settings();
    sets._exclusion_offsets = {0.5, 0.5};

    isobar pwave = solver.add_isobar<P_wave>(Delta, 1, id::Deck, "Deck", sets);
    
    // -----------------------------------------------------------------------
    timer timer;
    plotter plotter;

    timer.start();

    double smin = 0, smax = 2.0;
    double A = kinematics->A();
    double B = kinematics->B();
    double C = kinematics->C();
    double D = kinematics->D();

    plot p1 = plotter.new_plot();
    p1.set_legend(0.55, 0.6);
    p1.set_curve_points(1000);
    p1.set_xrange({smin, smax});
    p1.set_labels("#sigma   [GeV^{2}]", "#it{F}_{#Delta} (t, #it{m}_{3#pi}^{2} #; #sigma + #it{i}#epsilon)");
    p1.add_horizontal(0);
    p1.shade_region({A,C});
    p1.add_header("#minus #it{t}  = 0.1 GeV^{2}");
    p1.add_curve( {smin, smax}, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, "#Delta #times #Omega");
    p1.add_dashed({smin, smax}, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });
   
    std::array<std::string,4> labels = {"1st", "2nd", "3rd", "4th"};
    for (int i = 1; i <= N; i++)
    {
        solver.iterate();
        p1.add_curve( {smin, smax}, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, labels[i-1]);
        p1.add_dashed({smin, smax}, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); });
    };
    p1.save("deck_iterations.pdf");

    timer.stop();
    timer.print_elapsed();
};