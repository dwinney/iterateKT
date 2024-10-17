// The decay amplitude is decomposed into terms of one-variable functions
// These are given by the isobar class below
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
    inline isobar new_isobar(kinematics kin)
    {
        auto x = std::make_shared<A>(kin);
        return std::static_pointer_cast<raw_isobar>(x);
    };

    class raw_isobar
    {
        // -----------------------------------------------------------------------
        public:

        // Default constructor
        raw_isobar(kinematics xkin)
        : _kin(xkin)
        {};

        // -----------------------------------------------------------------------
        // Mandatory virtual methods which need to be overriden

        // Each isobar should have an identifying int (suggest implementing this with enums)
        virtual unsigned int id() = 0;

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
        complex omnes(double x, complex ieps = IEPS);

        // Evaluate the actual amplitude piece which usually involves a dispersion integral
        virtual inline complex evaluate(complex x){ return 0.; };

        // -----------------------------------------------------------------------
        // Utilities

        // Thing related to the options
        inline uint option(){ return _option; };
        virtual inline void set_option(uint x){ _option = x; };

 
        // -----------------------------------------------------------------------
        private:

        // Kinematics instance
        kinematics _kin;

        // Isobars name
        std::string _name = "isobar";

        // Simple id 
        unsigned int _id = 0;
        unsigned int _option = 0;

        // -----------------------------------------------------------------------
        protected:

        // Integrator settings
        int _integrator_depth = 5;

        // Interpolation settings
        double _interp_energy_low  = 5;    // interpolate from sth to this value
        int _interp_points_low     = 200;  // interpolate the above interval with this many points
        double _interp_energy_high = 1000; // then from _interp_energy_low to _interp_energy_high 
        int _interp_points_high    = 200;  // with this many points
    };

}; // namespace iterateKT

#endif // ISOBAR_HPP