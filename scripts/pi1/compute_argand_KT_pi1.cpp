// Export bare KT basis functions F_k(#sigma+i#epsilon) and bare Omnes #Omega(#sigma+i#epsilon)
// using the same setup as scripts/pi1/plot_isobars.cpp (P-wave pi1, four production mechanisms).
//
// No kinematic prefactor on F_k — rho_{\pi\pi}^{P-wave} = rho_CM * p_cm^2 and normalization are applied in plot_argand_KT_pi1.cpp.
//
// Output (under analysis/pi1/argand_KT_pi1/): one file per #it{m}_{3#pi} bin:
//   columns: sigma_GeV2 ReF0 ImF0 ... ReF3 ImF3 Omega_Re Omega_Im
//
// ------------------------------------------------------------------------------

#include "amplitude.hpp"
#include "utilities.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "settings.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

void compute_argand_KT_pi1()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    uint N_iter           = 9;
    uint asymptotic_power = 3;

    auto constant = [&](complex sigma) { return 1.; };
    auto gen_deck = [](double m3pi)
    {
        double t = -0.12;
        return [t, m3pi](complex sigma) { return pi1::deck(t, m3pi * m3pi, sigma); };
    };
    auto gen_bubble = [](double m3pi, double lam)
    {
        return [m3pi, lam](complex sigma) { return pi1::bubble(m3pi * m3pi, sigma, norm(lam)); };
    };

    timer timer;
    timer.start();

    std::vector<amplitude> amps;
    amps.emplace_back(new_amplitude<pi1>(new_kinematics(1.3, M_PION)));
    amps.emplace_back(new_amplitude<pi1>(new_kinematics(1.6, M_PION)));
    amps.emplace_back(new_amplitude<pi1>(new_kinematics(1.9, M_PION)));
    amps.emplace_back(new_amplitude<pi1>(new_kinematics(2.2, M_PION)));

    for (auto amp : amps)
    {
        double m3pi = amp->get_kinematics()->M();
        std::string name = "#it{m}_{3#pi} = " + to_string(m3pi) + " GeV";
        isobar pwave = amp->add_isobar<P_wave>(
            {constant, gen_bubble(m3pi, 0.05), gen_bubble(m3pi, 0.77), gen_deck(m3pi)},
            asymptotic_power, id::P_wave, name);
        amp->iterate(N_iter);
        timer.lap(name + " iterated!");
    }

    std::string out_dir = main_dir() + "/analysis/pi1/argand_KT_pi1/";
    std::filesystem::create_directories(out_dir);

    const std::array<double, 2> bounds = {0., 2.0};
    const int Npts                   = 1000;

    for (auto amp : amps)
    {
        double m3pi = amp->get_kinematics()->M();
        int tag     = static_cast<int>(std::round(m3pi * 1000.));

        isobar pwave = amp->get_isobar(id::P_wave);

        std::string path = out_dir + "m3pi_" + std::to_string(tag) + ".dat";
        std::ofstream out(path);
        out << std::scientific << std::setprecision(16);
        out << "# iterateKT bare F_k(sigma+i eps) and bare Omnes(sigma+i eps); rho and norm in plot script\n";
        out << "# matches plot_isobars.cpp (P-wave pi1)\n";
        out << "# m3pi_GeV " << m3pi << "\n";
        out << "# sigma_min_GeV2 " << bounds[0] << " sigma_max_GeV2 " << bounds[1] << " N " << Npts << "\n";
        out << "# columns: sigma_GeV2 ReF0 ImF0 ReF1 ImF1 ReF2 ImF2 ReF3 ImF3 Omega_Re Omega_Im\n";

        for (int i = 0; i < Npts; i++)
        {
            double sigma = bounds[0] + (bounds[1] - bounds[0]) * double(i) / double(Npts - 1);
            out << sigma;
            for (int k = 0; k < 4; k++)
            {
                complex F = pwave->basis_function(k, sigma + IEPS);
                out << " " << std::real(F) << " " << std::imag(F);
            }
            complex Oz = pwave->omnes(sigma + IEPS);
            out << " " << std::real(Oz) << " " << std::imag(Oz) << "\n";
        }
    }

    timer.stop();
    timer.print_elapsed();
}
