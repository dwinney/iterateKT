// The decay amplitude is decomposed into terms of one-variable functions
// These are given by the isobar class below. These will contain collections of
// iterations which contain the solutions at each step
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef ISOBAR_HPP
#define ISOBAR_HPP

#include <memory>
#include "utilities.hpp"
#include "kinematics.hpp"
#include "iteration.hpp"
#include "settings.hpp"
#include "basis.hpp"
#include "timer.hpp"
#include <Math/Interpolator.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/gauss.hpp>

namespace iterateKT
{
    // Forward declare for the typedef below
    class raw_isobar;

    // Define isobars only as pointers
    using isobar  = std::shared_ptr<raw_isobar>;

    // This function serves as our "constructor"
    template<class A>
    inline isobar new_isobar(kinematics kin, subtractions subs, uint maxsub, settings sets = settings() )
    {
        auto x = std::make_shared<A>(kin, subs, maxsub, sets);
        return std::static_pointer_cast<raw_isobar>(x);
    };

    class raw_isobar
    {
        // -----------------------------------------------------------------------
        public:

        // Default constructor
        raw_isobar(kinematics xkin, subtractions subs, uint maxsub, settings sets) 
        : _kinematics(xkin), _settings(sets), _subtractions(subs)
        { 
            // When we have "unsubtracted" we assume we do have one but no polynomial
            _max_sub = (maxsub == 0) ? 1 : maxsub;
            initialize();
        };

        // Each isobar should have an identifying int (suggest implementing this with enums)
        //  and a string name for human readable id
        inline void set_name(std::string x){ _name = x; };
        inline std::string name(){ return _name; }; 

        // -----------------------------------------------------------------------
        // Mandatory virtual methods which need to be overriden

        virtual id get_id() = 0;

        // The power (p q)^n that appears in angular momentum barrier factor
        // This determines the type of matching required at pseudothreshold
        virtual unsigned int singularity_power() = 0;

        // Elastic phase shift which provides the intial guess
        virtual double phase_shift(double s) = 0;

        // Kernel which appears in the angular average
        virtual complex ksf_kernel(id iso_id, complex s, complex t) = 0;

        // ----------------------------------------------------------------------- 
        // Things related to dispersion integrals and such

        // Evaluate the Omnes function (on the first sheet) for the given phaseshift. 
        // This is the homogenous solution and the start of our iterative procedure
        complex omnes(complex x);

        // Output a basis_function from a given iteration
        // Without an iter_id we just take the latest iteration
        complex basis_function(unsigned int iter_id, unsigned int basis_id, complex x);
        complex basis_function(unsigned int basis_id, complex x);

        // Take in an array of isobars and use their current state to calculate the next disc
        basis_grid calculate_next(std::vector<isobar> & previous_list);
        
        // Use the specified kernel to calculate the angualar average
        // over the pinocchio path
        complex pinocchio_integral(unsigned int basis_id, double s, std::vector<isobar> & previous_list);

        // Combine all the basis functions with their parameters
        inline complex evaluate(complex s)
        {
            complex sum = 0;
            for (int i = 0; i < _subtractions->N_basis(); i++) 
            sum += _subtractions->get_par(i)*basis_function(i, s);
            return sum;
        };

        // -----------------------------------------------------------------------
        // Utilities

        // Thing related to the options
        inline         uint option()          { return _option; };
        virtual inline void set_option(uint x){ _option = x; };

        // Flag used for internal debugging
        inline void set_debug(uint x){ _debug = x; };

        // Grab a pointer to a specific iteration of an isobar
        iteration get_iteration(uint id){ return _iterations[id];};

        // -----------------------------------------------------------------------
        protected:

        // Kinematics instance
        kinematics _kinematics;

        // Integrator and interpolator settings
        settings _settings;
        
        // Given a phase shift, this is the LHC piece (numerator)
        inline double LHC(double s)
        { 
            if (!_lhc_interpolated) interpolate_lhc();
            return (s > _kinematics->sth()) ? _lhc.Eval(s) : 0.; 
        };
        
        // Calculate the angular integral along a straight line
        // Bounds arguments should be {t_minus, t_plus, ieps perscription}
        complex linear_segment(unsigned int basis_id, std::array<double,3> bounds, double s, std::vector<isobar> & previous_list);
        // Calculate the integral along the curved secment of pinocchio's head
        complex curved_segment(unsigned int basis_id, double s, std::vector<isobar> & previous_list);

        // Save interpolation of the discontinuity calculated elsewhere into the list of iterations
        inline void save_iteration(basis_grid & grid){ _iterations.push_back(new_iteration( _kinematics, grid, _settings)); };

        // -----------------------------------------------------------------------
        private:

        // Private method only accessible to raw_amplitude
        friend class raw_amplitude;

        // Overal option flag 
        uint _option = 0;

        // IDs
        std::string _name = "isobar";

        // Debugging flag
        uint _debug = 0;
        bool debug(uint x){ return (_debug == x); };

        // Saved vector of iterations
        std::vector<iteration> _iterations;

        // Number of subtractions
        unsigned int _max_sub = 1, _n_basis = 1;
        subtractions _subtractions;

        // Initialize the 'zeroth' iteration by evaluating just the omnes function
        void initialize();
        std::vector<double> _s_list, _s_around_pth;
        complex _ieps;

        // Not properly initialized until lhc interpolation is saved
        void interpolate_lhc();
        ROOT::Math::Interpolator _lhc;
        bool _lhc_interpolated = false;
    };

}; // namespace iterateKT

#endif // ISOBAR_HPP