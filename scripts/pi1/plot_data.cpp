// Visualize the (complex) data from COMPASS collaboration [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2108.01744
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "plotter.hpp"
#include "fitter.hpp"

#include "COMPASS/data.hpp"

void plot_data()
{
    using namespace iterateKT;

    // -----------------------------------------------------------------------
    // Operating options

    // Data file
    std::string data_file = "dalitz_m3pi_bin_number_22_tBin_3.json";
    
    // Bounds for the axes
    std::array<double,2> xy_bounds = {0, 1.7}, z_bounds = {-1200, 1200};
    // Labels for the axes
    std::string xlabel = "#sigma_{b} [GeV^{2}]", ylabel =  "#sigma_{c} [GeV^{2}]";

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    auto  reim_data = COMPASS::parse_JSON_ReIm(data_file);
    auto  abs_data  = COMPASS::parse_JSON(data_file);
    
    // Need a kinematics instance to draw the dalitz boundary
    kinematics kin = new_kinematics(abs_data._extras["m3pi"], M_PION);

    // -----------------------------------------------------------------------
    // Plot results
    
    plotter plotter;
    
    // Next we make plots the raw data
    plot2D rep = kin->new_dalitz_plot(plotter);
    rep.set_Nbins(reim_data[0]._extras["Nbins"]);
    rep.set_data(reim_data[0]);
    rep.set_title("Real Part");
    rep.set_labels(xlabel, ylabel);
    rep.set_ranges(xy_bounds, xy_bounds);

    plot2D imp = kin->new_dalitz_plot(plotter);
    imp.set_Nbins(reim_data[1]._extras["Nbins"]);
    imp.set_data(reim_data[1]);
    imp.set_title("Imaginary Part");
    imp.set_labels(xlabel, ylabel);
    imp.set_ranges(xy_bounds, xy_bounds);

    plot2D abp = kin->new_dalitz_plot(plotter);
    abp.set_Nbins(abs_data._extras["Nbins"]);
    abp.set_data(abs_data);
    abp.set_title("Absolute Value");
    abp.set_labels(xlabel, ylabel);
    abp.set_ranges(xy_bounds, xy_bounds, {-EPS, 530});
      
    // Combine them all in one file
    plotter.combine({3,1}, {rep, imp, abp}, "data.pdf");
};