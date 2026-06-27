// Take in any number of m3pi and t bins (using COMPASS indexing) then iterate
// and export the KT solutions for each. 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Instituto de Ciencias Nucleares (ICN)
//               Universidad Nacional Autonoma de Mexico (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://inspirehep.net/literature/1898933
// ------------------------------------------------------------------------------


#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"
#include "COMPASS_pi1/fitter.hpp"

void bulk_iterate()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    // -----------------------------------------------------------------------
    // Operating options

    // How many times to iterate
    int niterate = 10;

    // Range of bins to consider (following COMPASS numbering) 
    int min = 11, max = 49;

    // Where to put files
    std::string export_path  = main_dir() + "/analysis/COMPASS_pi1/basis_functions/";

    // Prefix to label output files with
    std::string file_prefix  = "CD"; /* contact & Deck */

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution  
    
    
    // Grab sub-array of bin values
    std::vector<double> m_bins;
    for (int i = min; i <= max; i++) m_bins.push_back(COMPASS::m_bins[i-11]);
    
    // Set up our amplitude 
    amplitude amp  = new_amplitude<pi1_binned>(nullptr, std::make_tuple(m_bins, COMPASS::t_bins, niterate));
    
    // Check that our file directory exists
    std::filesystem::create_directory(export_path);
    // and export
    amp->export_solution(export_path + file_prefix); 
};