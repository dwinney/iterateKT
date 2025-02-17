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

        // -----------------------------------------------------------------------
        // Isobar management
        
        // Retrieve an isobar with index i
        inline isobar get_isobar(id i)
        { 
            for (auto f : _isobars) if (i == f->get_id()) return f;
            return error("amplitude::get_isobar", "Index out of scope!", nullptr);
        };

        // With just a single int nsub, we assume we have all subtraction coefficients to order nsub-1
        template<class T>
        inline void add_isobar(uint nsub, id id, std::string name = "isobar", settings sets = T::default_settings())
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

        // Else you can pass a vector of uints with the orders which get 
        template<class T>
        inline void add_isobar(std::vector<uint> poly, uint nsub, id id, std::string name = "isobar", settings sets = T::default_settings())
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
                    return;
                };
            };
            _isobars.push_back(new_iso);

            // Need to update _subtractions to know about new polynomials
            for (auto order : poly)
            {
                _subtractions->_ids.push_back(new_iso->get_id());
                _subtractions->_powers.push_back(order);
            };
        };

        // Access the full vector of isobar pointers
        inline std::vector<isobar> get_isobars(){ return _isobars; };

        // Print to file necessary info to reconstruct isobars later
        void export_solution(std::string prefix);
      
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