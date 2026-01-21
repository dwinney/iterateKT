// Take in any number of m3pi and t bins (using COMPASS indexing) then iterate
// and export KT solutions for each
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

    int niterate = 10;
    std::string export_path  = "/scripts/pi1/basis_functions/";
    std::string file_prefix  = "CD"; /* Contact & Deck */

    // -----------------------------------------------------------------------
    // Check that our file directory exists

    std::filesystem::create_directory(main_dir()+export_path);

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution

    std::vector<double> m_bins = {1.3, 1.4, 1.5};
    // Set up our amplitude 
    amplitude amp  = new_amplitude<pi1_binned>(nullptr, std::make_tuple(m_bins, COMPASS::t_bins, niterate));
    amp->export_solution(export_path + file_prefix); 
};