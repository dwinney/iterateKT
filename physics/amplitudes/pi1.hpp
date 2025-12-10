// Amplitudes relevant for the decay of meson with JP = 1-+ into 3pi
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef PI1_AMPLITUDES_HPP
#define PI1_AMPLITUDES_HPP

#include <tuple>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "amplitude.hpp"
#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "timer.hpp"
#include"isobars/pi1.hpp"

namespace iterateKT
{ 
    inline settings default_settings()
    {
        settings sets;
        sets._exclusion_points        = 10;
        sets._infinitesimal           = 1E-8;
        sets._intermediate_energy     = 5;
        sets._cutoff                  = 20;
        sets._interpolation_offset    = 0.1;
        sets._interpolation_points    = {200, 10, 100};
        double xi_sth = 1E-3,  eps_sth = 1E-3;
        double xi_pth = 3E-3,  eps_pth = 4E-3;
        double xi_rth = 3E-1,  eps_rth = 3E-1;

        sets._exclusion_offsets   = {5E-2, 5E-2};
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};

        phase_args iso_1   = {"madrid/delta_11.dat", 1.69,  1, 1, 2};
        sets._phase_shifts = { {id::P_wave, iso_1}, {id::Contact, iso_1}, {id::Deck, iso_1} };

        return sets;
    };

    class pi1 : public raw_amplitude
    {
        public: 
        
        // Constructor
        pi1(kinematics kin) : raw_amplitude(kin) {};
        
        // Spin 1 decay so (2j+1) = 3
        inline double combinatorial_factor(){ return 3; };

        static constexpr double _mu  = 0.13957000;
        static constexpr double _mu2 = _mu*_mu;

        // ------------------------------------------------------------------------------
        // Things related to the inclusion of the bubble

        // P-wave projection of 2 particle phase space with BW factor
        // M2   -> total 3pi invariant mass
        // s    -> 2pi subsystem mass
        // lam2 -> cutoff mass squared (mass of first ignored particle)
        static inline complex bubble(complex M2, complex s, double lam2)
        {
            complex mu2 = complex(_mu2);
            complex q   = csqrt(kallen(M2, s, mu2))/2/csqrt(M2);
            complex rho = 2*q/csqrt(M2);

            // careful if we cross above the three-body cut
            // multiply by -1 to not change sign and stay on the same sheet
            bool above_3bcut = (real(s) >= real(M2)+_mu2);
            if  (above_3bcut){ q *= -1; rho *= -1; };

            complex z2 = q*q / lam2;
            complex blatt_weisskopf = 1./(1+z2);
            return z2*blatt_weisskopf*rho;
        };

        // ------------------------------------------------------------------------------
        // Things related to the inclusion of the Deck triangle

        static inline complex tau(complex t, complex M2, complex s, double z)
        {
            complex mu2 = complex(_mu2);
            complex p   = csqrt(kallen(M2, t, mu2))/2/csqrt(M2);
            complex q   = csqrt(kallen(M2, s, mu2))/2/csqrt(M2);
            complex rho = 2*q/csqrt(M2);

            // careful if we cross above the three-body cut
            // multiply by -1 to not change sign and stay on the same sheet
            bool above_3bcut = (real(s) >= real(M2)+_mu2);
            if  (above_3bcut){ q *= -1; rho *= -1; };
            return  mu2 + s - (M2+mu2-t)*(M2+s-mu2)/2/M2 + 2*p*q*z;
        };

        // 3P1 projection of the vanilla OPE
        // t    -> momentum transfer of (external) Pomeron
        // M2   -> total 3pi invariant mass
        // s    -> 2pi subsystem imvariant mass
        // mex2 -> exchanged particle mass
        static inline complex deck(complex t, complex M2, complex s, double mex2)
        {
            // Masses and momenta
            complex mu2 = complex(_mu2);
            complex p   = csqrt(kallen(M2, t, mu2))/2/csqrt(M2);
            complex q   = csqrt(kallen(M2, s, mu2))/2/csqrt(M2);
            complex rho = 2*q/csqrt(M2);

            // careful if we cross above the three-body cut
            // multiply by -1 to not change sign and stay on the same sheet
            bool above_3bcut = (real(s) >= real(M2)+_mu2);
            if  (above_3bcut){ q *= -1; rho *= -1; };
            // Momentum tranfer at costheta = 0
            complex t0  = tau(t, M2, s, 0); 
            // Angular argument
            complex z   = (mex2 - t0)/2/p/q;

            // Continution depends on the ieps used for s
            // and not that of z
            bool above_thr  = real(s) >= 4*_mu2;

            // Legendre of 2nd kind
            complex Q0;
            if (above_thr) Q0 = log(-(z+1)/(z-1))/2-I*sign(imag(s))*PI/2;
            else           Q0 = log( (z+1)/(z-1))/2;

            // Final discontinuity
            return rho*q/p*((1-z*z)*Q0+z);
        };

        // without specififying an exchange mass, assume we mean pion
        static inline complex deck(complex t, complex M2, complex s)
        { return deck(t, M2, s, _mu2); };     
        // adding a form factor is just the pion exchange minus the cutoff      
        static inline complex deck_with_FF(complex t, complex M2, complex s, double lam2)
        { return deck(t, M2, s, _mu2) - deck(t, M2, s, lam2); };

        // Assuming a pi- pi- pi+ decay and only P-waves
        // s = (pi- + pi+)^2 
        // t = (pi- + pi+)^2
        // u = (pi- + pi-)^2
        inline complex prefactor_s(id iso_id, complex s, complex t, complex u){ return csqrt(_kinematics->kibble(s,t,u)); };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u){ return - prefactor_s(iso_id, t, s, u); };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u){ return 0.; };
    };

    // ------------------------------------------------------------------------------
    // The following are container amplitudes which holds multiple copies of the above
    // but at different bins in production t for simultaneous fits

    enum class option : unsigned int { set_m3pibin, set_tbin, set_m3pibin_COMPASS };

    class pi1_across_tbins : public raw_amplitude
    {
        public:

        pi1_across_tbins(kinematics xkin, std::tuple<std::array<double,4>,uint> args)
        : raw_amplitude(xkin), _tvals(std::get<0>(args)), _niter(std::get<1>(args))
        {
            for (int i = 0; i < 4; i++) 
            { 
                _tbins.emplace_back(new_amplitude<pi1>(xkin));
                _tbins[i]->set_name("tbin " + to_string(i));
            };
            initialize();
        };

        inline void set_option(option opt, double x)
        {
            switch (opt)
            {
                case option::set_tbin:  _current = _tbins[int(std::round(x))]; break;
                default: return;
            };
        };
        inline uint N_pars(){ return _current->N_pars(); };
        inline void set_parameters(std::vector<complex> x){ _current->set_parameters(x); };
        inline complex evaluate(complex s, complex t, complex u){ return _current->evaluate(s, t, u); };

        private:

        // how many times to iterate by default
        int _niter = 10;

        // Save each of the 4 bins as a pointer
        std::array<double,4>   _tvals;
        std::vector<amplitude> _tbins;

        // This is the one that gets called
        amplitude _current;

        inline void initialize()
        {
            auto constant = [ ](complex sigma){ return complex(1.); };
            auto deck     = [&](double t)
            { 
                double m3pi2 = _kinematics->M2();
                return [m3pi2,t](complex sigma){ return pi1::deck(t, m3pi2, sigma);}; 
            };

            // Set up all the amplitudes
            for (int i = 0; i < 4; i++)
            {
                _tbins[i]->add_isobar<P_wave>({constant, deck(_tvals[i])}, 3, id::P_wave, "P-wave");
                _tbins[i]->iterate(_niter);
            };

            // At end place the first one in _current
            set_option(option::set_tbin, 0);
        };
    };

    // ------------------------------------------------------------------------------
    // Increasing complexity, this class allows simultaneous 2D fits to a bunch of
    // bins in both m3pi and t

    class pi1_binned : public raw_amplitude
    {
        public:

        // We require an overall xkin but this will be ignored since we
        // create individual instances for each needed amplitude
        pi1_binned(kinematics xkin, std::tuple<std::vector<double>,std::array<double,4>,uint> args)
        : raw_amplitude(xkin)
        {
            timer timer;
            timer.start();
            auto m3pi_vals = std::get<0>(args);
            auto t_vals    = std::get<1>(args);
            auto nint      = std::get<2>(args);
            for (auto m3pi : m3pi_vals) 
            { 
                _kins.emplace_back(new_kinematics(m3pi, M_PION));
                _m3pibins.emplace_back(new_amplitude<pi1_across_tbins>(_kins.back(), std::make_tuple(t_vals, nint)));
                timer.lap("initialized amplitude with m3pi = "+to_string(m3pi));
            };
            timer.stop(); 
            timer.print_elapsed();
            set_option(option::set_m3pibin, 0);
            set_option(option::set_tbin,    0);
        };

        // Most utilities should just pipe to whatever _current is pointed to 
        inline kinematics get_kinematics(){ return _current->get_kinematics(); };
        inline void set_parameters(std::vector<complex> x){ _current->set_parameters(x); };
        inline complex evaluate(complex s, complex t, complex u){ return _current->evaluate(s, t, u); };

        // Except the total number of pars which are cumulative
        inline uint N_pars(){ return _current->N_pars()*_m3pibins.size(); };

        // Maneuver which subamplitude we're pointing to
        inline void set_option(option opt, double x)
        {
            int ix = int(std::round(x));
            switch (opt)
            {
                case option::set_m3pibin_COMPASS:  
                {
                    _current = _m3pibins[find_COMPASS_bin(ix)]; 
                    _current->set_option(option::set_tbin, _current_tbin);
                    break;
                };
                case option::set_m3pibin:  
                {
                    _current = _m3pibins[ix]; 
                    _current->set_option(option::set_tbin, _current_tbin);
                    break;
                };
                case option::set_tbin:   
                {
                    _current->set_option(option::set_tbin, x); 
                    _current_tbin = ix;
                    break;
                };
                default: return;
            };
        };

        private: 

        // Store each m3pi requires its own kinematics and amplitude
        std::vector<kinematics> _kins;
        std::vector<amplitude>  _m3pibins;

        // The way we navigate bins, we want to keep track of which tbin we're looking at
        uint _current_tbin = 0;

        // pointer to the current amplitude to be evaluated
        amplitude _current; 

        // Map the COMPASS bin codes [11 - 48] to the indexes of our saved amplitudes
        inline int find_COMPASS_bin(uint bin)
        {
            // Bins are 40 MeV wide and start at 0.96 GeV
            double M = 0.96 + (bin-11)*0.04;
            for (int i = 0; i < _m3pibins.size(); i++)
            {
                double Mi = _kins[i]->M();
                if (are_equal(Mi, M)) return i;
            };
            return -1;
        };
    };
}; // namespace iterateKT 

#endif // PI1_AMPLITUDES_HPP