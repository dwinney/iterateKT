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
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace iterateKT
{
    // We need a way to differentiate different channels
    // since isobars describe physics in one channel
    // at a time
    enum channel{S_CHANNEL, T_CHANNEL, U_CHANNEL};

    // Forward declare for the typedef below
    class raw_isobar;

    // Define isobars only as pointers
    using isobar  = std::shared_ptr<raw_isobar>;

    // This function serves as our "constructor"
    template<class A>
    inline isobar new_isobar(kinematics kin, int nsub, settings sets = settings() )
    {
        auto x = std::make_shared<A>(kin, nsub, sets);
        return std::static_pointer_cast<raw_isobar>(x);
    };

    class raw_isobar
    {
        // -----------------------------------------------------------------------
        public:

        // Default constructor
        raw_isobar(kinematics xkin, int nsub, settings sets) : _kinematics(xkin), _settings(sets)
        { 
            set_max_subtraction(nsub); 
            initialize();
        };

        // -----------------------------------------------------------------------
        // Mandatory virtual methods which need to be overriden

        // Each isobar should have an identifying int (suggest implementing this with enums)
        virtual unsigned int id()  = 0;
        virtual std::string name() = 0; // as well as a string name for human readable id

        // Elastic phase shift which provides the intial guess
        virtual double phase_shift(double s) = 0;

        // In the full amplitude, user must define the prefactor associated with
        // each channel. This is where symmetries, kinematic factors, and angular polynomials,
        // are implemented;
        virtual complex prefactor(channel chan, complex s, complex t, complex u) = 0;

        // ----------------------------------------------------------------------- 
        // Things related to dispersion integrals and such

        // Evaluate the Omnes function (on the first sheet) for the given phaseshift. 
        // This is the homogenous solution and the start of our iterative procedure
        complex omnes(complex x);

        // Output a basis_function from a given iteration
        // Without an iter_id we just take the latest iteration
        complex basis_function(unsigned int iter_id, unsigned int basis_id, complex x);
        complex basis_function(unsigned int basis_id, complex x);

        // Evaluate the full isobar combining the basis functions and coefficients
        // Specify an iter_id or just eval the latest one
        complex evaluate(unsigned int iter_id, complex s);
        complex evaluate(complex s);

        // -----------------------------------------------------------------------
        // Utilities

        // Thing related to the options
        inline uint option(){ return _option; };
        virtual inline void set_option(uint x){ _option = x; };

        // -----------------------------------------------------------------------
        protected:

        friend class amplitude;

        // Take in an array of isobars and use their current state to calculate the next disc
        std::array<std::vector<double>,3> calculate_next(std::vector<isobar> previous);

        // Save interpolation of the discontinuity calculated elsewhere into the list of iterations
        inline void save_iteration(std::array<std::vector<double>,3> disc)
        {
            _iterations.push_back(new_iteration(_max_sub, disc, _settings));
        }

        // -----------------------------------------------------------------------
        private:

        // Kinematics instance
        kinematics _kinematics;

        // Integrator and interpolator settings
        settings _settings;

        // Overal option flag 
        int _option = 0;

        // Saved vector of iterations
        std::vector<iteration> _iterations;

        // Number of subtractions
        unsigned int _max_sub = 1;
        void set_max_subtraction(int n)
        {
            std::string message = "Isobar initiated with 0 subtraction will be ignored and initialized with 1 instead.";
            _max_sub = (n == 0) ? error(message, 1) : n;

            // Initialize each subtraction coefficient to 0
            for (int i = 0; i < n; i++) _subtraction_coeffs.push_back(0.); 
        };

        // Initialize the 'zeroth' iteration by evaluating just the omnes function
        inline void initialize(){ _iterations.push_back(new_iteration(_max_sub)); };

        // Subtraction coefficients
        std::vector<complex> _subtraction_coeffs;

    };

}; // namespace iterateKT

#endif // ISOBAR_HPP