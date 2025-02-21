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
    inline amplitude new_amplitude(kinematics kin, std::string id = "amplitude")
    {
        auto x = std::make_shared<A>(kin, id);
        return std::static_pointer_cast<raw_amplitude>(x);
    };


    class raw_amplitude : public solver
    {
        // -----------------------------------------------------------------------
        public:

        // Define only the masses here. 
        // The amplitude structure from quantum numbers will come later
        raw_amplitude(kinematics xkin, std::string id) : solver(xkin)
        {};
   
        // Evaluate the full amplitude.
        virtual complex evaluate(complex s, complex t, complex u);

        // Factor to divide by in width calculation
        virtual double  helicity_factor(){ return 1; };

        // Need to specify how to combine the isobars into the full amplitude
        virtual complex prefactor_s(id iso_id, complex s, complex t, complex u){ return 0.; };
        virtual complex prefactor_t(id iso_id, complex s, complex t, complex u){ return 0.; };
        virtual complex prefactor_u(id iso_id, complex s, complex t, complex u){ return 0.; };

        // Calculate widths in the physical decay region
        double differential_width(double s, double t);
        double differential_width(double s);
        double width();

        // -----------------------------------------------------------------------
        // Utilities

        // Set and get the string id 
        inline void set_name(std::string name){ _name = name; };
        inline std::string name(){ return _name; };

        // Pass an option flag and so something. By default we dont do anything 
        // Can be overloaded to do whatever you want 
        virtual inline void set_option(option opt){ return; }; 

        inline void set_parameters( std::vector<complex> pars)
        {
            if (pars.size() != _subtractions->N_basis())
            {
                warning("set_parameters", "Parameter vector of unexpected size!");
                return;
            }
            _subtractions->_values = pars;
        };

        // Number of free parameters (used by fitters)
        inline uint N_pars(){ return _subtractions->N_basis(); };

        inline std::vector<complex> get_pars() { return _subtractions->_values; };

        // -----------------------------------------------------------------------
        private:

        // Id string to identify the amplitude with
        std::string _name = "amplitude";

        // (inverse of) prefactors for the differential width
        inline double prefactors(){ return 32*pow(2*PI*_kinematics->M(),3); }
    };
}; // namespace iterateOKT

#endif // AMPLITUDE_HPP