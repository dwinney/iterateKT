// Minimal structure to solve the KT equations
// This should be used if only isobars are important and not any amplitude
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <memory>
#include "kinematics.hpp"
#include "utilities.hpp"
#include "basis.hpp"
#include "isobar.hpp"

namespace iterateKT
{
    // Forward declare function to be specified at runtime
    settings default_settings();

    class solver 
    {
        // -----------------------------------------------------------------------
        public: 

        solver(kinematics xkin)
        : _kinematics(xkin), _subtractions(std::make_shared<raw_subtractions>())
        {};

        // Calculate one iteration of the KT equations
        void iterate();
        inline void iterate(unsigned int N){ for (int i = 0; i < N; i++) iterate(); };

        // Iterate N times with nice little terminal messages
        void timed_iterate(unsigned int N);

        // -----------------------------------------------------------------------
        // Isobar management
        
        // Retrieve an isobar with index i
        inline isobar get_isobar(id i)
        { 
            for (auto f : _isobars) if (i == f->get_id()) return f;
            return error("amplitude::get_isobar", "Index out of scope!", nullptr);
        };

        // Grab an isobar and in the same step give it name.
        // Chances are you already named it in the output variable so just pass that name to the isobar itself too
        inline isobar get_isobar(id i, std::string name)
        { 
            isobar x = get_isobar(i); 
            if(x != nullptr) x->set_name(name);
            return x;            
        };

        // Search for isobar through its name
        inline isobar get_isobar(std::string name)
        { 
            for (auto isobar : _isobars)
            {
                if (isobar->name() == name) return isobar;
            };
            return error("amplitude::get_isobar", "Cannot find isobar named + "+name+ "!", nullptr);
        };

        // Else you can pass a vector of uints with the orders which get 
        template<class T>
        inline isobar add_isobar(std::vector<std::function<complex(complex)>> driving_terms, uint nsub, 
                               id id, std::string name = "isobar", settings sets = default_settings())
        { 
            isobar_args args;
            args._kin      = _kinematics;
            args._id       = id;
            args._subs     = _subtractions;
            args._name     = name;
            args._sets     = sets;
            args._maxsub   = nsub;
            isobar new_iso = new_isobar<T>(args);

            // Check for uniqueness
            for (auto old_iso : _isobars)
            {
                if (old_iso->get_id() == new_iso->get_id())
                {
                    warning("add_isobar", "Attempted to add an isobar with non-unique identifier!");
                    return nullptr;
                };
            };
            _isobars.push_back(new_iso);

            // Need to update _subtractions to know about new polynomials
            for (auto & dts : driving_terms)
            {
                _subtractions->_ids.push_back(new_iso->get_id());
                _subtractions->_drivers.push_back(std::move(dts));
            };

            return new_iso;
        };

        template<class T>
        inline isobar add_isobar(std::function<complex(complex)> driving_term, uint nsub, id id, 
                                                                                        std::string name = "isobar", 
                                                                                        settings sets = default_settings())
        { 
            return add_isobar<T>(std::vector<std::function<complex(complex)>>{driving_term}, nsub, id, name, sets);
        };

        // Else you can pass a vector of uints with the orders which get 
        template<class T>
        inline isobar add_isobar(std::vector<uint> poly, uint nsub, id id, std::string name = "isobar", 
                                                                         settings sets = default_settings())
        { 
            std::vector<std::function<complex(complex)>> driving_terms;
            for (auto power : poly)
            {
                if (power == 0) driving_terms.push_back([power](complex x){ return 1.; });
                else driving_terms.push_back([power](complex x){ return pow(x,power); });
            };
            return add_isobar<T>(driving_terms, nsub, id, name, sets);
        };
    
        // With just a single int nsub, we assume we have all subtraction coefficients to order nsub-1
        template<class T>
        inline isobar add_isobar(uint nsub, id id, std::string name = "isobar", 
                                                 settings sets = default_settings())
        { 
            std::vector<uint> pows;
            for (int i = 0; i < nsub; i++) pows.push_back(i);
            return add_isobar<T>(pows, nsub, id, name, sets);
        };

        // Access the full vector of isobar pointers
        inline std::vector<isobar> get_isobars(){ return _isobars; };

        // Print to file necessary info to reconstruct isobars later
        void export_solution(std::string prefix, uint precision = 12);

        // output the saved kinematics
        inline kinematics get_kinematics(){ return _kinematics; };
      
        // -----------------------------------------------------------------------
        protected: 
        
        // Kinematics object, contains all masses, angles, etc
        kinematics _kinematics;

        // Store isobars here to be called later
        std::vector<isobar> _isobars;

        // Store info regarding the subtraction polynomials
        subtractions _subtractions;
    };

}; // namespace iterateKT

#endif // SOLVER_HPP