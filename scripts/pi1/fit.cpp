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
#include "COMPASS/data.hpp"

void fit()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // -----------------------------------------------------------------------
    // Operating options

    int bin_number = 15;  // which m3pi bin to fit
    int Niter      = 5;   // Number of KT iterations
    int Nsub       = 2;   // Number of subtractions

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution

    // Import our data set first so we can know the m3pi bin
    data_set data   = COMPASS::parse_JSON("dalitz_m3pi_bin_number_"+to_string(bin_number)+"_tBin_3.json");
    double m3pi     = data._extras["m3pi"];

    // Set up general kinematics so everything knows masses
    kinematics kin = new_kinematics(m3pi, M_PION);
    
    // Set up our amplitude 
    amplitude amp = new_amplitude<pi1>(kin, "π₁ → 3π");

    // We only have one isobar, which we import here
    amp->add_isobar<P_wave>(Nsub, id::P_wave);

    // Iterate Niter times
    amp->timed_iterate(Niter);

    // -----------------------------------------------------------------------
    // Set up fitter

    // These vectors should be same size as Nsub above
    std::vector<std::string> all_labels = {"alpha", "beta", "gamma"};
    std::vector<std::string> par_labels(all_labels.begin(), all_labels.begin() + Nsub);
    std::vector<complex> initial_guess(Nsub, 1.0);

    // Add data
    fitter<COMPASS::fit> fitter(amp);
    fitter.add_data(data);
    
    fitter.set_parameter_labels(par_labels);
    fitter.fix_argument("alpha", 0.); // Fix overall phase 

    fitter.do_fit(initial_guess);

    // -----------------------------------------------------------------------
    // Plot results

    plotter plotter;

    std::array<double,2> bounds = {0, kin->pth()+0.1};
    std::string xlabel = "#sigma_{b}  [GeV^{2}]", ylabel =  "#sigma_{c}  [GeV^{2}]";

    // Plot the amplitude
    plot2D p1 = amp->plot_dalitz(plotter);
    p1.set_palette(kBird);
    p1.set_labels(xlabel, ylabel);
    p1.set_ranges(bounds, bounds);
    
    // Finally calculatet the chi2 per bin
    std::vector<double> pull;
    for (int i = 0; i < data._N; i++)
    {
        double s1 = data._x[i], s2 = data._y[i];
        complex model = amp->evaluate(s1, s2);

        double fcn = (is_zero(data._dz[i])) ? 0. : (std::abs(model) - data._z[i]) / data._dz[i];
        pull.push_back(fcn);
    };
    double max_pull = *std::max_element(pull.begin(), pull.end());

    plot2D p2 = kin->new_dalitz_plot(plotter);
    p2.set_palette(kTemperatureMap);
    p2.set_data({data._x, data._y, pull});
    p2.set_labels(xlabel, ylabel);
    p2.set_ranges(bounds, bounds, {-max_pull, max_pull});

    // Combine them all in one file
    plotter.combine({2,1}, {p1,p2}, "fit_results.pdf");
    
    std::vector<double> bins, ends, model_in_bin; 
    double max_z    = *std::max_element(data._z.begin(), data._z.end());
    for (int i = 0; i < data._N; i++) 
    {
        bins.push_back(i);
        
        double s = data._x[i], t = data._y[i];
        complex M = amp->evaluate(s, t);
        model_in_bin.push_back( abs(M) );
    };
    
    double n = data._N / 12;
    std::vector<plot>   bin_plots;
    for (int i = 0; i < 12; i++)
    {
        plot p = plotter.new_plot();
        p.add_data(bins, {data._z, data._dz});
        p.add_curve(bins, model_in_bin);
        p.set_ranges({n*i, n*(i+1)}, {0, max_z});
        bin_plots.push_back(p);
    };

    plotter.combine({4,3}, bin_plots, "bins.pdf");
};