// Fit KT amplitudes for π1(1600) decay with only P-wave in [1]
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
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "plotter.hpp"
#include "fitter.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"
#include "COMPASS/data.hpp"

void fit()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // -----------------------------------------------------------------------
    // Operating options

    // Data file
    std::string data_file = "dalitz_m3pi_bin_number_22_tBin_3.json";

    int Niter = 5; // Number of KT iterations
    int Nsub  = 1; // Number of subtractions
    
    // These vectors should be same size as Nsub above
    std::vector<std::string> par_labels = {"a", "b"};
    std::vector<complex> initial_guess  = { 1.,  1.};
    
    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    // Import our data set first so we can know the m3pi bin
    data_set data   = COMPASS::parse_JSON(data_file);
    double m3pi     = data._extras["m3pi"];

    // Set up general kinematics so everything knows masses
    kinematics kin = new_kinematics(m3pi, M_PION);
    
    // Set up our amplitude 
    amplitude amp = new_amplitude<pi1>(kin, "π₁ → 3π");

    // We only have one isobar, which we import here
    amp->add_isobar<P_wave>(Nsub, id::P_wave);

    // Iterate
    divider();
    print("Solving KT with " + to_string(Nsub) + " subtractions and " + to_string(Niter) + " iterations:");
    line();
    timer timer; 
    timer.start();
    for (int i = 0; i < Niter; i++)
    {
        amp->iterate();
        timer.lap("iteration " + to_string(i+1));
    };
    timer.stop();
    line();

    timer.print_elapsed();
    divider(); 

    // -----------------------------------------------------------------------
    // Set up fitter

    fitter<COMPASS::fit> fitter(amp);
    fitter.set_parameter_labels(par_labels);
    fitter.add_data(data);
    fitter.do_fit(initial_guess);

    // -----------------------------------------------------------------------
    // Plot results

    plotter plotter;

    std::array<double,2> bounds = {0, 1.7};
    std::string xlabel = "#sigma_{1} [GeV^{2}]", ylabel =  "#sigma_{2} [GeV^{2}]";

    // Plot the amplitude
    plot2D p1 = amp->plot_dalitz(plotter, "[GeV^{2}]");
    p1.set_title("|Fit|");
    p1.set_ranges(bounds, bounds, {-EPS, 1000});
    
    // Finally calculatet the chi2 per bin
    std::vector<double> chi2;
    for (int i = 0; i < data._N; i++)
    {
        double s1 = data._x[i], s2 = data._y[i];
        complex model = amp->evaluate(s1, s2);

        double fcn = (is_zero(data._dz[i])) ? 0. : std::norm(data._z[i] - std::abs(model))/data._dz[i];
        chi2.push_back(fcn/fitter.dof());
    };

    plot2D p2 = kin->new_dalitz_plot(plotter);
    p2.set_palette(kTemperatureMap);
    p2.set_data({data._x, data._y, chi2});
    p2.set_title("#chi^{2} / dof");
    p2.set_labels(xlabel, ylabel);
    p2.set_ranges(bounds, bounds, {-EPS, 7});

    // Combine them all in one file
    plotter.combine({2,1}, {p1,p2}, "results.pdf");
};