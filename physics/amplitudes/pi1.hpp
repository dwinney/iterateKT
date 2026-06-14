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
        sets._infinitesimal           = 1E-7;
        sets._intermediate_energy     = 2.0;
        sets._cutoff                  = 20;
        sets._interpolation_offset    = 0.001;
        sets._interpolation_points    = {200, 8, 100};
        double xi_sth = 1E-3,  eps_sth = 1E-3;
        double xi_pth = 1E-3,  eps_pth = 4E-2;
        double xi_rth = 2E-1,  eps_rth = 6E-2;

        sets._exclusion_offsets   = {7E-2, 2E-1};
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};

        std::string path = main_dir() + "/physics/phase_shifts/";
        phase_args iso_1   = {path+"madrid/delta_11.dat", 1.69,  1, 1, 2};
        sets._phase_shifts = { {id::P_wave, iso_1}, {id::Contact, iso_1}, {id::Deck, iso_1} };

        return sets;
    };

    class pi1 : public raw_amplitude
    {
        public: 
        
        // Constructor
        pi1(kinematics kin) : raw_amplitude(kin) {};
        
        // Spin 1 decay so (2j+1) = 3
        inline double combinatorial_factor(){ return 1.; };

        static constexpr double _mu  = 0.13957000;
        static constexpr double _mu2 = _mu*_mu;

        // ------------------------------------------------------------------------------
        // Things related to the inclusion of the Deck triangle

        // 3P1 projection of the vanilla OPE
        // t     -> momentum transfer of (external) Pomeron
        // m3pi2 -> total 3pi invariant mass 
        // sigma -> 2pi subsystem imvariant mass
        static inline complex deck(double t, double m3pi2, complex sigma)
        {
            double eps = imag(sigma);
            bool almost_real = abs(eps) < 1E-5;

            // Thresholds and other special points
            double sth  =  norm(2*_mu);
            double pth  =  norm( sqrt(m3pi2) - _mu );
            double rth  =  norm( sqrt(m3pi2) + _mu );
            double tbc  =  m3pi2 + _mu2;
            double pole = (m3pi2 - _mu2)*(m3pi2 - t + _mu2)/(m3pi2 + t - _mu2);
            double exchange_mass = _mu2;

            // Momenta
            complex p   = csqrt(kallen( m3pi2, t,     _mu2))/2/csqrt(m3pi2);
            complex q   = csqrt(kallen( m3pi2, sigma, _mu2))/2/csqrt(m3pi2);
            // careful if we cross above the m3pi2 + mu2 we have a vertical cut
            bool above_tbc = real(sigma) >= tbc;
            if  (above_tbc) q *= -1;
            
            // Momentum tranfer at costheta = 0
            complex t0  = _mu2 + sigma - (m3pi2+_mu2-t)*(m3pi2+sigma-_mu2)/2/m3pi2; 
            // Angular argument
            complex z   = (exchange_mass - t0)/2/p/q;

            // Legendre of the second kind
            complex Q0;
            Q0 = (imag(z) > 0) ? log((z+1)/(z-1))/2 : (log((1+z)/(1-z))-I*PI)/2;
            if (eps < 0) Q0 += I*PI;

            // Final result with spin factors
            return ((1-z*z)*Q0+z)/p/q;
        };    

        // When we want the integrated width to just be the integrated intenstiy 
        // and so we remove the flux prefactors to match COMPASS definition
        inline double prefactors(){ return 1; };

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

    enum class option : unsigned int { set_mbin, set_tbin, set_mbin_COMPASS };

    class pi1_across_tbins : public raw_amplitude
    {
        public:

        // With no extra int assume we dont iterate
        pi1_across_tbins(kinematics xkin, std::tuple<std::array<double,4>> args)
        : raw_amplitude(xkin), _tvals(std::get<0>(args))
        {
            for (int i = 0; i < 4; i++) 
            { 
                _tbins.emplace_back(new_amplitude<pi1>(xkin));
                _tbins[i]->set_name("tbin_" + to_string(i));
            };
            initialize(0);
        };

        pi1_across_tbins(kinematics xkin, std::tuple<std::array<double,4>,uint> args)
        : raw_amplitude(xkin), _tvals(std::get<0>(args))
        {
            for (int i = 0; i < 4; i++) 
            { 
                _tbins.emplace_back(new_amplitude<pi1>(xkin));
                _tbins[i]->set_name("tbin_" + to_string(i));
            };
            initialize(std::get<1>(args));
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
        inline complex evaluate_in_dalitz(double s, double t){ return _current->evaluate_in_dalitz(s, t); };

        // Observables
        inline double width(){ return _current->width() ; };

        // Export solution iterates over the four tbins exporting each one
        inline void export_solution(std::string path, uint precision)
        {          
            std::string file = main_dir() + path;
            for (int i = 0; i < 4; i++) _tbins[i]->export_solution(path, precision);
        };

        // Similar with import
        inline void import_solution(std::string path)
        {
            for (int i = 0; i < 4; i++)
            {
                std::string file = path + "_t_" + to_string(-_tvals[i]) + ".dat";
                _tbins[i]->get_isobars()[0]->import_iteration<3>(file);
            };
        };

        inline void precompute_dalitz(uint N)
        {
            for (auto bin : _tbins) bin->precompute_dalitz(N);
        };  

        private:

        // Save each of the 4 bins as a pointer
        std::array<double,4>   _tvals;
        std::vector<amplitude> _tbins;

        // This is the one that gets called
        amplitude _current;

        inline void initialize(uint niter)
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
                _tbins[i]->add_isobar<P_wave>({constant, constant, deck(_tvals[i])}, 3, id::P_wave, "t_"+to_string(-_tvals[i]));
                _tbins[i]->iterate(niter);
            };

            // At end place the first one in _current
            set_option(option::set_tbin, 0);
        };
    };

    // ------------------------------------------------------------------------------
    // Increasing complexity, this class allows simultaneous 2D fits to a bunch of bins in both m3pi and t

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
            auto m_vals  = std::get<0>(args);
            auto t_vals  = std::get<1>(args);
            auto nint    = std::get<2>(args);
            for (auto m3pi : m_vals) 
            { 
                _mbins.emplace_back(new_amplitude<pi1_across_tbins>(new_kinematics(m3pi, M_PION), std::make_tuple(t_vals, nint)));
                timer.lap("initialized amplitude with m3pi = "+to_string(m3pi));
            };
            timer.stop(); 
            timer.print_elapsed();
            set_option(option::set_mbin, 0);
            set_option(option::set_tbin, 0);
        };

        // With no last uint, we assume uniterated
        pi1_binned(kinematics xkin, std::tuple<std::vector<double>,std::array<double,4>> args)
        : raw_amplitude(xkin)
        {
            auto m_vals  = std::get<0>(args);
            auto t_vals  = std::get<1>(args);
            for (auto m3pi : m_vals) 
            { 
                _mbins.emplace_back(new_amplitude<pi1_across_tbins>(new_kinematics(m3pi, M_PION), std::make_tuple(t_vals)));
            };
            set_option(option::set_mbin, 0);
            set_option(option::set_tbin, 0);
        };

        // Most utilities should just pipe to whatever _current is pointed to 
        inline kinematics get_kinematics(){ return _current->get_kinematics(); };
        inline void set_parameters(std::vector<complex> x){ _current->set_parameters(x); };
        inline complex evaluate(complex s, complex t, complex u){ return _current->evaluate(s, t, u); };
        inline complex evaluate_in_dalitz(double s, double t){ return _current->evaluate_in_dalitz(s, t); };

        // Observables
        inline double width(){ return _current->width() ; };

        // Except the total number of pars which are cumulative
        inline uint N_pars(){ return _current->N_pars()*_mbins.size(); };

        // Maneuver which subamplitude we're pointing to
        inline void set_option(option opt, double x)
        {
            int ix = int(std::round(x));
            switch (opt)
            {
                case option::set_mbin_COMPASS:  
                {
                    _current = _mbins[find_COMPASS_bin(ix)]; 
                    _current->set_option(option::set_tbin, _current_tbin);
                    break;
                };
                case option::set_mbin:  
                {
                    _current = _mbins[ix]; 
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

        // Export solution now simply iterates over mbins making sure we name things correctly
        inline void export_solution(std::string path, uint precision)
        {
            for (auto bin : _mbins)
            {
                double m = bin->get_kinematics()->M();
                std::string file = path + "_m_" + to_string(m);
                bin->export_solution(file, precision);
            }
        };

        // Similar with import
        inline void import_solution(std::string path)
        {
            for (auto bin : _mbins)
            {
                double m = bin->get_kinematics()->M();
                std::string new_path = path + "_m_" + to_string(m);
                bin->import_solution(new_path);
            };
        };

        inline void precompute_dalitz(uint N)
        {
            for (auto bin : _mbins) bin->precompute_dalitz(N);
            return;
        };

        protected: 

        // Store each m3pi requires its own kinematics and amplitude
        std::vector<amplitude>  _mbins;

        // The way we navigate bins, we want to keep track of which tbin we're looking at
        uint _current_tbin = 0;

        // pointer to the current amplitude to be evaluated
        amplitude _current; 

        // Map the COMPASS bin codes [11 - 48] to the indexes of our saved amplitudes
        inline int find_COMPASS_bin(uint bin)
        {
            // Bins are 40 MeV wide and start at 0.96 GeV
            double M = 0.96 + (bin-11)*0.04;
            for (int i = 0; i < _mbins.size(); i++)
            {
                double Mi = _mbins[i]->get_kinematics()->M();
                if (are_equal(Mi, M)) return i;
            };
            return -1;
        };
    };
}; // namespace iterateKT 

#endif // PI1_AMPLITUDES_HPP