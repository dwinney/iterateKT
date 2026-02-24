// Fit KT amplitudes for π1(1600) decay with only P-wave to data from [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://inspirehep.net/literature/1898933
// ------------------------------------------------------------------------------

#include <algorithm>
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "plotter.hpp"
#include "fitter.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"
#include "COMPASS_pi1/fitter.hpp"
#include "COMPASS_pi1/data.hpp"

void single_fit()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // -----------------------------------------------------------------------
    // Operating options

    int m3pibin    = 22;  // which m3pi bin to fit 22, 27 & 32
    int tbin       = 0;   // which t bin to fit
    int Niter      = 10;  // Number of KT iterations

    // Import our data set first so we can know the m3pi bin
    data_set data   = COMPASS::parse_JSON(m3pibin, tbin);
    double m3pi     = data._extras["m3pi"];
    double t        = data._extras["t"];

    // Contact piece gets just constant as driving term
    auto   constant = [&](complex sigma){return 1.;};
    auto   linear   = [&](complex sigma){return sigma;};
    auto   quad     = [&](complex sigma){return sigma*sigma;};
    auto   bubble   = [&](complex sigma){return pi1::bubble(m3pi*m3pi, sigma, norm(0.2));};
    auto   deck     = [&](complex sigma){return pi1::deck(t, m3pi*m3pi, sigma);};

    std::vector<std::function<complex(complex)>> driving_terms = {constant, deck};

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution

    // Set up general kinematics so everything knows masses
    kinematics kin = new_kinematics(m3pi, M_PION);
    
    // Set up our amplitude 
    amplitude amp  = new_amplitude<pi1>(kin);
    amp->set_name("π₁ → 3π");

    // Add isobar using the above function as our driving term
    isobar pwave   =  amp->add_isobar<P_wave>(driving_terms,  3, id::P_wave, "Deck");

    // Iterate Niter times
    amp->timed_iterate(Niter);

    // -----------------------------------------------------------------------
    // Set up fitter

    // These vectors should be same size as Nsub above
    std::vector<complex> initial_guess;
    // for (auto x : driving_terms) initial_guess.push_back(1.0);
    // initial_guess = {1009.56026045, complex(-1677.79215508,180.617166393)};
    initial_guess = {319.8, 755.8*exp(I*3.13)};

    // Add data
    fitter<COMPASS::fit_single_bin> fitter(amp, "Combined");
    fitter.set_tolerance(0.00001E3);
    fitter.set_print_level(3);
    fitter.add_data(data);
    
    fitter.make_real("par[0]"); 
    fitter.do_fit(initial_guess);

    // -----------------------------------------------------------------------
    // Plot results

    plotter plotter;

    std::array<double,2> bounds = {0, kin->pth()+0.1};
    std::string xlabel = "#sigma_{#it{a}}  [GeV^{2}]", ylabel =  "#sigma_{#it{b}}  [GeV^{2}]";

    // Plot the amplitude
    plot2D p1 = amp->plot_dalitz(plotter);
    p1.set_palette(kBird);
    p1.set_labels(xlabel, ylabel);
    p1.set_ranges(bounds, bounds);
    p1.save("dalitz.pdf");
    
    // Finally calculatet the chi2 per bin
    std::vector<double> pull, bin_i;
    for (int i = 0; i < data._N; i++)
    {
        bin_i.push_back(i);
        double s1 = data._x[i], s2 = data._y[i];
        complex model = amp->evaluate(s1, s2, amp->get_kinematics()->Sigma() - s1 - s2);

        double fcn = (is_zero(data._dz[i])) ? 0. : (std::abs(model) - data._z[i]) / data._dz[i];
        pull.push_back(fcn);
    };
    double max_pull = *std::max_element(pull.begin(), pull.end());

    plot2D p2 = kin->new_dalitz_plot(plotter);
    p2.set_Nbins(data._extras["Nbins"]);
    p2.set_palette(kTemperatureMap);
    p2.set_data({data._x, data._y, pull});
    p2.set_labels(xlabel, ylabel);
    p2.set_ranges(bounds, bounds, {-max_pull, max_pull});
    p2.save("pull_2D.pdf");

    plot p3 = plotter.new_plot();
    p3.add_data(bin_i, pull, dot(jpacColor::DarkGrey));
    p3.add_horizontal(0);
    p3.set_labels("Bin number", "Pull");
    p3.set_ranges( {0, bin_i.back()}, {-6,6});
    p3.save("pull_1D.pdf");
};