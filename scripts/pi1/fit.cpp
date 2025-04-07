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
    auto data   = COMPASS::parse_JSON(data_file);

    double m3pi = data[0]._extras["m3pi"]; // 3pi decay mass in GeV
    double mt   = data[1]._extras["t"];    // Production -t

    int Niter = 0; // Number of KT iterations
    int Nsub  = 1; // Number of subtractions

    // These vectors should be same size as Nsub above
    std::vector<std::string> par_labels = {"a", "b"};
    std::vector<complex> initial_guess  = { 1.,  1.};

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution

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
    fitter.add_data<2>(data);
    fitter.do_fit(initial_guess);

    // -----------------------------------------------------------------------
    // Plot results

    plotter plotter;

    // ------------------------------------------------------------
    // 3D PLOTS 

    // These two plots are of the amplitude
    auto plots = amp->make_plots(plotter, "[GeV^{2}]");
    plots[0].set_title("Re(Fit)");
    plots[0].set_ranges( {0, 1.7}, {0, 1.7}, {-1200, 1200});
    plots[1].set_title("Im(Fit)"); 
    plots[1].set_ranges( {0, 1.7}, {0, 1.7}, {-1200, 1200});
    
    // Next we make plots the raw data
    plot2D rep = kin->new_dalitz_plot(plotter);
    rep.set_data(data[0]);
    rep.set_title("Re(Data)");
    rep.set_labels("#sigma_{1} [GeV^{2}]", "#sigma_{2} [GeV^{2}]");
    rep.set_ranges( {0, 1.7}, {0, 1.7}, {-1200, 1200});

    plot2D imp = kin->new_dalitz_plot(plotter);
    imp.set_data(data[1]);
    imp.set_title("Im(Data)");
    imp.set_labels("#sigma_{1} [GeV^{2}]", "#sigma_{2} [GeV^{2}]");
    imp.set_ranges( {0, 1.7}, {0, 1.7}, {-1200, 1200});

    // Finally we can calculate and plot the difference
    // between the two
    std::vector<double> dre, dim;
    for (int i = 0; i < data[0]._N; i++)
    {
        double s1 = data[0]._x[i], s2 = data[0]._y[i];
        complex model = amp->evaluate(s1, s2);
        dre.push_back(data[0]._z[i] - real(model));
        dim.push_back(data[1]._z[i] - imag(model));
    };

    plot2D drep = kin->new_dalitz_plot(plotter);
    drep.set_palette(kTemperatureMap);
    drep.set_data({data[0]._x, data[0]._y, dre});
    drep.set_title("Re(Data - #it{A})");
    drep.set_labels("#sigma_{1} [GeV^{2}]", "#sigma_{2} [GeV^{2}]");
    drep.set_ranges( {0, 1.7}, {0, 1.7});

    plot2D dimp = kin->new_dalitz_plot(plotter);
    dimp.set_palette(kTemperatureMap);
    dimp.set_data({data[0]._x, data[0]._y, dim});
    dimp.set_title("Im(Data - #it{A})");
    dimp.set_labels("#sigma_{1} [GeV^{2}]", "#sigma_{2} [GeV^{2}]");
    dimp.set_ranges( {0, 1.7}, {0, 1.7});

    // Combine them all in one file
    plotter.combine({2,3}, {rep, imp, plots[0], plots[1], drep, dimp}, "dalitz.pdf");
};