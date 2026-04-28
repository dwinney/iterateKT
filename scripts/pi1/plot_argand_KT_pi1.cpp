// Argand diagrams from compute_argand_KT_pi1.cpp (bare F_k and Omega).
//
// Applies rho_{\pi\pi}(\sigma)\, p_{\mathrm{cm}}^{2}(\sigma) (P-wave barrier times CM rho),
// then normalizes each KT trajectory to +i at sigma = M_RHO^2. Omnes normalized separately.
//
// Reads: analysis/pi1/argand_KT_pi1/m3pi_*.dat
//
// ------------------------------------------------------------------------------

#include "constants.hpp"
#include "utilities.hpp"
#include "plotter.hpp"
#include "colors.hpp"
#include "data_set.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void plot_argand_KT_pi1()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // rho_{pipi}^{P-wave}(sigma) = rho_CM(sigma) * p_cm^2(sigma); rho_CM = sqrt(1 - sth/sigma),
    // p_cm^2 = (sigma - sth)/4 for equal-mass pi pi (sth = 4 mpi^2).
    auto rho_P_wave_pi_pi = [](double sigma)
    {
        const double sth = 4. * M_PION * M_PION;
        if (sigma <= sth)
            return 0.;
        double rho_cm = std::sqrt(1. - sth / sigma);
        double p_sq = (sigma - sth) / 4.;
        return rho_cm * p_sq;
    };

    std::array<int, 4> tags = {1300, 1600, 1900, 2200};
    std::array<double, 4> m3pi_vals = {1.3, 1.6, 1.9, 2.2};
    std::array<std::string, 4> basis_labels = {
        "Contact",
        "Bubble (Lambda = 200 MeV)",
        "Bubble (Lambda = 770 MeV)",
        "Deck (t = -0.12 GeV^{2})",
    };

    std::array<jpacColor, 4> basis_colors = {
        jpacColor::Blue, jpacColor::Red, jpacColor::Green, jpacColor::Orange};

    std::string dir = main_dir() + "/analysis/pi1/argand_KT_pi1/";

    plotter plotter_obj;
    std::vector<plot> panels;

    const std::array<double, 2> xrange = {-0.7, 0.7};
    const std::array<double, 2> yrange = {-0.2, 1.2};
    const double s_rho                    = M_RHO * M_RHO;

    for (size_t b = 0; b < tags.size(); b++)
    {
        std::string data_path = dir + "m3pi_" + std::to_string(tags[b]) + ".dat";
        auto cols             = import_data<11>(data_path);
        if (cols[0].empty())
        {
            std::cerr << "Missing precomputed data: " << data_path
                      << "\nRun compute_argand_KT_pi1.cpp first.\n";
            return;
        }

        // Plot only up to sigma_plot_max (GeV^2). High-sigma tail is numerically noisy for some bases
        // (e.g. Bubble Lambda = 200 MeV). Optional extra trim of last points after the cutoff.
        const double sigma_plot_max = 1.42;
        const size_t trim_extra     = 40;

        size_t n_plot = 0;
        for (size_t row = 0; row < cols[0].size(); ++row)
        {
            if (cols[0][row] <= sigma_plot_max)
                n_plot = row + 1;
            else
                break;
        }
        if (n_plot > trim_extra + 5)
            n_plot -= trim_extra;
        if (n_plot < 2)
        {
            n_plot = cols[0].size() < 50 ? cols[0].size() : 50;
        }

        size_t ir = 0;
        for (size_t row = 1; row < cols[0].size(); row++)
        {
            if (std::fabs(cols[0][row] - s_rho) < std::fabs(cols[0][ir] - s_rho))
                ir = row;
        }

        const double rho_r = rho_P_wave_pi_pi(cols[0][ir]);
        complex norms[4];
        for (int k = 0; k < 4; k++)
        {
            complex F_at =
                complex(cols[static_cast<size_t>(1 + 2 * k)][ir], cols[static_cast<size_t>(2 + 2 * k)][ir]);
            complex q_at = rho_r * F_at;
            if (std::abs(q_at) < 1e-14)
                std::cerr << "Warning: rho_{pipi}^{P-wave} * F nearly zero at rho row for basis " << k << ", m3pi="
                          << m3pi_vals[b] << "\n";
            norms[k] = I / q_at;
        }

        complex Om_at =
            complex(cols[9][ir], cols[10][ir]);
        if (std::abs(Om_at) < 1e-14)
            std::cerr << "Warning: Omnes nearly zero at rho row, m3pi=" << m3pi_vals[b] << "\n";
        complex norm_O = I / Om_at;

        plot p = plotter_obj.new_plot();

        std::ostringstream title;
        title << "m_{3#pi} = " << std::fixed << std::setprecision(2) << m3pi_vals[b] << " GeV";
        p.set_plot_title(title.str());
        p.set_ranges(xrange, yrange);
        // ROOT TLatex: displays rho_{pipi}^P as requested (P-wave pipi kinematic factor).
        p.set_labels(
            "Re[F_{k} #rho_{#pi#pi}^{P}]",
            "Im[F_{k} #rho_{#pi#pi}^{P}]");
        p.set_legend(0.35, 0.38);
        p.set_legend_spacing(0.036);

        p.add_horizontal(0.);
        p.add_vertical(0.);

        std::vector<double> oxre, oyim;
        oxre.reserve(n_plot);
        oyim.reserve(n_plot);
        for (size_t row = 0; row < n_plot; row++)
        {
            complex Oz = norm_O * complex(cols[9][row], cols[10][row]);
            oxre.push_back(std::real(Oz));
            oyim.push_back(std::imag(Oz));
        }
        p.add_curve(oxre, oyim, solid(jpacColor::DarkGrey, "Omnes (Madrid)"));

        for (int k = 0; k < 4; k++)
        {
            std::vector<double> xre, yim;
            xre.reserve(n_plot);
            yim.reserve(n_plot);
            for (size_t row = 0; row < n_plot; row++)
            {
                double rhoP = rho_P_wave_pi_pi(cols[0][row]);
                complex F =
                    complex(cols[static_cast<size_t>(1 + 2 * k)][row], cols[static_cast<size_t>(2 + 2 * k)][row]);
                complex z = norms[k] * rhoP * F;
                xre.push_back(std::real(z));
                yim.push_back(std::imag(z));
            }
            entry_style style = solid(basis_colors[k], basis_labels[k]);
            p.add_curve(xre, yim, style);
        }

        std::vector<double> px = {0.};
        std::vector<double> py = {1.};
        p.add_data(px, py, dot(jpacColor::DarkGrey, ""));

        panels.push_back(p);
    }

    plotter_obj.combine({2, 2}, panels, "argand_KT_pi1.pdf");
}
