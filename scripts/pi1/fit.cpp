// Twice subtracted KT amplitudes for π1(1600) decay with only P-wave in [1]
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

    int Niter = 5; // Number of KT iterations
    int Nsub  = 2; // Number of subtractions

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
    amp->iterate(Niter);
    
    // -----------------------------------------------------------------------
    // Set up fitter

    fitter<COMPASS::fit> fitter(amp);
    fitter.set_parameter_labels(par_labels);
    fitter.add_data<2>(data);
    fitter.do_fit(initial_guess);

    // -----------------------------------------------------------------------
    // Plot results

    plotter plotter;
    auto plots = amp->make_plots(plotter, "[GeV^{2}]");
    
    // Plot Dalitz Plots

    plot2D rep = kin->new_dalitz_plot(plotter);
    rep.set_data(data[0]);
    rep.set_title("Re(Data)");
    rep.set_labels("#sigma_{1} [GeV^{2}]", "#sigma_{2} [GeV^{2}]");

    plot2D imp = kin->new_dalitz_plot(plotter);
    imp.set_data(data[1]);
    imp.set_title("Im(Data)");
    imp.set_labels("#sigma_{1} [GeV^{2}]", "#sigma_{2} [GeV^{2}]");

    plotter.combine({2,2}, {rep, imp, plots[0], plots[1]}, "dalitz.pdf");

    // Plot the isobar compares to the Omness
    isobar pwave = amp->get_isobar(id::P_wave);
    complex norm = amp->get_parameters().front();

    plot p = plotter.new_plot();
    p.set_legend(0.7, 0.7);
    p.set_labels("#sigma  [GeV^{2}]", "#it{f}#kern[-0.4]{_{1}}#kern[-1]{^{1}}(#sigma)");
    p.add_curve({0, 1.5}, [&](double s){ return real(pwave->omnes(s+IEPS)); },          solid( jpacColor::Blue, "Omn#grave{e}s"));
    p.add_curve({0, 1.5}, [&](double s){ return imag(pwave->omnes(s+IEPS)); },          dotted(jpacColor::Blue));
    p.add_curve({0, 1.5}, [&](double s){ return real(pwave->evaluate(s+IEPS)/norm); },  solid( jpacColor::Red,  "KT ("+to_string(Nsub)+"-sub)"));
    p.add_curve({0, 1.5}, [&](double s){ return imag(pwave->evaluate(s+IEPS)/norm ); }, dotted(jpacColor::Red));
    p.save("isobar.pdf");
};