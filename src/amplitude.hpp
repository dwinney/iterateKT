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

namespace iterateKT
{
    // Forward declare for the typedef below
    class raw_amplitude;

    // Define isobars only as pointers
    using amplitude  = std::shared_ptr<raw_amplitude>;

    // This function serves as our "constructor"
    template<class A>
    inline amplitude new_amplitude(kinematics kin, std::string id = "amplitude")
    {
        auto x = std::make_shared<A>(kin, id);
        return std::static_pointer_cast<raw_amplitude>(x);
    };


    class raw_amplitude
    {
        // -----------------------------------------------------------------------
        public:

        // Define only the masses here. 
        // The amplitude structure from quantum numbers will come later
        raw_amplitude(kinematics xkin, std::string id)
        : _kinematics(xkin), _id(id), _subtractions(std::make_shared<raw_subtractions>())
        {};

        // Evaluate the full amplitude. This will 
        virtual complex evaluate(complex s, complex t, complex u);

        // Need to specify how to combine the isobars into the full amplitude
        virtual complex prefactor_s(id iso_id, complex s, complex t, complex u){ return 0.; };
        virtual complex prefactor_t(id iso_id, complex s, complex t, complex u){ return 0.; };
        virtual complex prefactor_u(id iso_id, complex s, complex t, complex u){ return 0.; };

        // Calculate one iteration of the KT equations
        void iterate();
        inline void iterate(unsigned int N){ for (int i = 0; i < N; i++) iterate(); };
        
        // -----------------------------------------------------------------------
        // Utilities

        // Retrieve or set the option flag.
        // set_option can be overloaded if you want to do more than just save it      
        inline int option(){ return _option; };
        virtual inline void set_option(int x){ _option = x; for (auto f : _isobars) f->set_option(x); };

        inline void set_parameters( std::vector<complex> pars)
        {
            if (pars.size() != _subtractions->N_basis())
            {
                warning("set_parameters", "Parameter vector of unexpected size!");
                return;
            }
            _subtractions->_values = pars;
        };

        // -----------------------------------------------------------------------
        // Isobar management
        
        // Retrieve an isobar with index i
        inline isobar get_isobar(id i)
        { 
            for (auto f : _isobars) if (i == f->get_id()) return f;
            return error("amplitude::get_isobar", "Index out of scope!", nullptr);
        };

        // Load up a new isobar
        template<class T>
        inline void add_isobar(int nsub, settings sets = T::default_settings())
        { 
            isobar new_iso = new_isobar<T>(_kinematics, _subtractions, nsub, sets);
            
            // Check for uniqueness
            for (auto old_iso : _isobars)
            {
                if (old_iso->get_id() == new_iso->get_id())
                {
                    warning("add_isobar", "Attempted to add an isobar with non-unique identifier!");
                    return;
                };
            };
            _isobars.push_back(new_iso);

            // Need to update _subtractions to know about new polynomials
            for (int i = 0; i < nsub; i++)
            {
                _subtractions->_ids.push_back(new_iso->get_id());
                _subtractions->_powers.push_back(i);
            };
        };

        // Access the full vector of isobar pointers
        inline std::vector<isobar> get_isobars(){return _isobars;};

        // -----------------------------------------------------------------------

        private:

        // Kinematics object, contains all masses, angles, etc
        kinematics _kinematics;

        // Store isobars here to be called later
        std::vector<isobar> _isobars;

        // Store info regarding the subtraction polynomials
        subtractions _subtractions;

        // Options flag
        int _option;

        // Id string to identify the amplitude with
        std::string _id = "amplitude";
    };
}; // namespace iterateOKT

#endif // AMPLITUDE_HPP