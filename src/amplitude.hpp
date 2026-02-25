// Top level class which defines the 1 -> 3 decay process. 
// The template feeds in all the relevant user info for a specific process.
// So far we require all equal mass particles
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef AMPLITUDE_HPP
#define AMPLITUDE_HPP

#include <memory>
#include "kinematics.hpp"
#include "utilities.hpp"
#include "basis.hpp"
#include "isobar.hpp"
#include "solver.hpp"
#include "plotter.hpp"
#include "plot2D.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace iterateKT
{
    // Forward declare the option type
    enum class option : unsigned int;

    // Forward declare for the typedef below
    class raw_amplitude;

    // Define isobars only as pointers
    using amplitude  = std::shared_ptr<raw_amplitude>;

    // This function serves as our "constructor"
    template<class A=raw_amplitude>
    inline amplitude new_amplitude(kinematics kin)
    {
        auto x = std::make_shared<A>(kin);
        return std::static_pointer_cast<raw_amplitude>(x);
    };

    // This function serves as our "constructor" with some extra needed arg
    template<class A=raw_amplitude, typename B>
    inline amplitude new_amplitude(kinematics kin, B extra_args)
    {
        auto x = std::make_shared<A>(kin, extra_args);
        return std::static_pointer_cast<raw_amplitude>(x);
    };

    class raw_amplitude : public solver
    {
        // -----------------------------------------------------------------------
        public:

        // Define only the masses here. 
        // The amplitude structure from quantum numbers will come later
        raw_amplitude(kinematics xkin) : solver(xkin) {};
   
        // Evaluate the full amplitude.
        virtual complex evaluate(complex s, complex t, complex u);

        // Factor to divide by in width calculation
        virtual double  combinatorial_factor(){ return 1; };

        // Need to specify how to combine the isobars into the full amplitude
        virtual complex prefactor_s(id iso_id, complex s, complex t, complex u){ return 0.; };
        virtual complex prefactor_t(id iso_id, complex s, complex t, complex u){ return 0.; };
        virtual complex prefactor_u(id iso_id, complex s, complex t, complex u){ return 0.; };

        // Calculate widths in the physical (symmetric) decay region
        virtual double differential_width(double s, double t);
        virtual double differential_width(double s);
        virtual double width();

        // (inverse of) prefactors for the differential width
        virtual inline double prefactors()
        { return 32*pow(2*PI*_kinematics->M(),3)*combinatorial_factor(); }

        // -----------------------------------------------------------------------
        // Utilities

        // Set and get the string id 
        inline void set_name(std::string name){ _name = name; };
        inline std::string name(){ return _name; };

        // Pass an option flag and so something. By default we dont do anything 
        // Can be overloaded to do whatever you want 
        virtual inline void set_option(option opt)          { return; }; 
        virtual inline void set_option(option opt, double x){ return; }; 

        virtual inline void set_parameters( std::vector<complex> pars)
        {
            if (pars.size() != _subtractions->N_basis())
            {
                warning("set_parameters", "Parameter vector of unexpected size!");
                return;
            }
            _subtractions->_values = pars;
        };

        // Number of free parameters (used by fitters)
        virtual inline uint N_pars(){ return _subtractions->N_basis(); };

        // Get a vector of the currently saved subtraction parameters
        virtual inline std::vector<complex> get_pars()       { return _subtractions->_values; };
        virtual inline std::vector<complex> get_parameters() { return get_pars(); };

        // Return a pointer to the interally saved kinematics instance
        virtual inline kinematics get_kinematics() { return _kinematics;}

        // Importing solutions from file are typically trickier than exporting
        // So we kep this as an overrideable function which does nothing by defualt
        // the importing procedure should be written by user
        virtual inline void import_solution(std::string path){ return; };

        // Some amplitudes contain pointers to other amplitudes this can be specified to allowed one to access
        // the currently pointed to subamplitude. By default this will return nullptr
        virtual inline amplitude get_current(){ return nullptr; };

        // -----------------------------------------------------------------------
        // Automate making plots of the amplitude

        // Make standard plots in the entire decay region
        // output[0] = Re(A)
        // output[1] = Im(A)
        std::vector<plot2D> plot_ReIm(plotter & pltr, std::string units = "", int N = 100);

        // Only plot |A|
        plot2D plot_dalitz(plotter & pltr, std::string units = "", int N = 100);

        // -----------------------------------------------------------------------
        // Calculate Daltiz plot parameters from amplitude
        // By default we use the conventions and notation of the PDG 
        // see ``Dalitz Plot Parameters for K -> 3pi decays" in RPP
        // Output in order {g, h, j, k, f}

        // The optional arguments are for the center of the dalitz plot and normalization:
        // X = (t - s )/m[0]
        // Y = (u - s0)/m[1]

        std::array<double,5> get_dalitz_parameters(double eps, double s0, std::array<double,2> ms);

        // Use the s0 from _kinematics
        std::array<double,5> get_dalitz_parameters(double eps = 1E-5, std::array<double,2> ms = {M_PION*M_PION, M_PION*M_PION})
        {
            return get_dalitz_parameters(eps, _kinematics->s0(), ms);
        };

        // -----------------------------------------------------------------------
        private:

        // Id string to identify the amplitude with
        std::string _name = "amplitude";
    };
}; // namespace iterateOKT

#endif // AMPLITUDE_HPP