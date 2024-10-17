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
#include "isobar.hpp"

namespace iteratedOKT
{
    class amplitude
    {
        // -----------------------------------------------------------------------
        public:

        // Define only the masses here. 
        // The amplitude structure from quantum numbers will come later
        amplitude(kinematics xkin)
        : _kin(xkin)
        {};

        // Evaluate the full amplitude. This will 
        complex operator()(complex s, complex t, complex u);

        // -----------------------------------------------------------------------
        // Utilities

        // Retrieve or set the option flag.
        // set_option can be overloaded if you want to do more than just save it      
        inline int option(){ return _option; };
        virtual inline void set_option(int x){ _option = x; for (auto f : _isobars) f->set_option(x); };

        // -----------------------------------------------------------------------
        // Isobar management
        
        // Retrieve an isobar with index i
        inline isobar get_isobar(unsigned int i)
        { 
            if (i > _isobars.size()) return error("amplitude::get_isobar", "Index out of scope!", nullptr);
            return _isobars[i]; 
        };

        // Load up a new isobar
        template<class T>
        inline void add_isobar(){ _isobars.push_back(new_isobar<T>(_kin)); };

        // -----------------------------------------------------------------------

        private:

        // Kinematics object, contains all masses, angles, etc
        kinematics _kin;

        // Store isobars here to be called later
        std::vector<isobar> _isobars;

        // Options flag
        int _option;

    };
}; // namespace iterateOKT

#endif // AMPLITUDE_HPP