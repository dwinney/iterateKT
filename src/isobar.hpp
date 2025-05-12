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
#include "phase_shift.hpp"
#include "basis.hpp"
#include "data_set.hpp"
#include "timer.hpp"
#include <Math/Interpolator.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace iterateKT
{
    // Forward declare for the typedef below
    class raw_isobar;

    // Define isobars only as pointers
    using isobar  = std::shared_ptr<raw_isobar>;

    // All the things you need to initialize an isobar can get lumped into a struct
    struct isobar_args
    {
        kinematics _kin;
        id _id;
        subtractions _subs;
        uint _maxsub;
        std::string _name;
        settings _sets;
    };

    // This function serves as our "constructor"
    template<class A>
    inline isobar new_isobar(isobar_args args)
    {
        return std::static_pointer_cast<raw_isobar>(std::make_shared<A>(args));
    };

    class raw_isobar
    {
        // -----------------------------------------------------------------------
        public:

        // Default constructor
        raw_isobar(isobar_args args) 
        : _kinematics(args._kin), _settings(args._sets), _subtractions(args._subs), 
                                  _id(args._id), _name(args._name)
        { 
            // When we have "unsubtracted" we assume we do have one but no polynomial
            _max_sub = (args._maxsub == 0) ? 1 : args._maxsub;
            initialize();
        };

        raw_isobar(isobar_args args, phase_args phase_args) 
        : _kinematics(args._kin), _settings(args._sets), _subtractions(args._subs), 
                                  _id(args._id), _name(args._name), _delta(phase_args)
        { 
            // When we have "unsubtracted" we assume we do have one but no polynomial
            _max_sub = (args._maxsub == 0) ? 1 : args._maxsub;
            initialize();
        };

        // -----------------------------------------------------------------------
        // Mandatory virtual methods which need to be overriden

        // The power (p q)^2*n+1 that appears in angular momentum barrier factor
        // This determines the type of matching required at pseudothreshold
        virtual unsigned int angular_momentum() = 0;

        // Elastic phase shift which provides the initial guess
        // By default we assume the _delta was imported with an interpolation
        virtual double phase_shift(double s){ return _delta(s); };

        // Kernel which appears in the angular average
        virtual complex ksf_kernel(id iso_id, complex s, complex t){ return 0.; };

        // This is an extra function to check if to calculate an inhomogeneity at all
        // Will return false unless overriden
        virtual bool calculate_inhomogeneity(){ return true; };

        // ----------------------------------------------------------------------- 
        // Things related to dispersion integrals and such

        // Evaluate the Omnes function (on the first sheet) for the given phaseshift. 
        // This is the homogenous solution and the start of our iterative procedure
        complex omnes(complex x);

        // Given a phase shift, this is the LHC piece (numerator)
        inline double LHC(double s)
        { 
            if (!_lhc_interpolated) interpolate_lhc();
            return (s > _kinematics->sth()) ? _lhc.Eval(s) : 0.; 
        };
                
        // Output a basis_function from a given iteration
        // Without an iter_id we just take the latest iteration
        complex basis_function(unsigned int iter_id, unsigned int basis_id, complex x);
        complex basis_function(unsigned int basis_id, complex x);

        // Calculate the first and second derivatives of basis function at the subraction point (s=0)
        template<uint order>
        inline complex basis_derivative(uint basis_id, double x, double eps)
        {
            auto f = [this,basis_id](double s){ return basis_function(basis_id, s); };
            return central_difference_derivative<complex>(order, f, x, eps);
        };

        // Take in an array of isobars and use their current state to calculate the next disc
        basis_grid calculate_next(std::vector<isobar> & previous_list);
        
        // Use the specified kernel to calculate the angualar average
        // over the pinocchio path
        complex pinocchio_integral(unsigned int basis_id, double s, std::vector<isobar> & previous_list);

        // Combine all the basis functions with their parameters
        inline complex evaluate(complex s)
        {
            complex sum = 0;
            for (int i = 0; i < _subtractions->N_basis(); i++) sum += _subtractions->get_par(i)*basis_function(i, s);
            return sum;
        };

        // -----------------------------------------------------------------------
        // Utilities
        
        // Import an iteration that was previously saved to file
        template<int N>
        void import_iteration(std::string filename)
        {
            auto imported = import_data<2*N+1>(filename, true);
            _s_list       = std::move(imported[0]);

            if (_s_list.back() < _settings._cutoff) 
            fatal("import_isobar", 
                  "Cutoff set to higher than imported data! Will cause interpolation errors.");

            basis_grid grid;
            grid._n_singularity = 2*angular_momentum()+1;
            grid._s_list        = _s_list;
            grid._s_around_pth  = _s_around_pth;
            
            for (int i = 0; i < N; i++)
            {
                grid._re_list.push_back(std::move(imported[1+2*i]));
                grid._im_list.push_back(std::move(imported[2+2*i]));
            };
            save_iteration(grid);
        };

        // Flag used for internal debugging
        inline void set_debug(uint x){ _debug = x; };
        // Grab pointer to latest iteration
        iteration get_iteration(){ return _iterations.back();};
        // Grab a pointer to a specific iteration of an isobar
        iteration get_iteration(uint id){ return _iterations[id];};

        // Each isobar should have an identifying int (suggest implementing this with enums)
        //  and a string name for human readable id
        inline void set_name(std::string x){ _name = x; };
        inline std::string name(){ return _name; }; 

        // Only way to publically access the isobars id
        inline id get_id(){ return _id; };

        // -----------------------------------------------------------------------
        protected:

        // Id of the isobar, see definition in basis.hpp
        id _id;

        // Kinematics instance
        kinematics _kinematics;

        // Integrator and interpolator settings
        settings _settings;

        // Calculate the angular integral along a straight line
        // Bounds arguments should be {t_minus, t_plus, ieps perscription}
        complex linear_segment(unsigned int basis_id, std::array<double,3> bounds, double s, std::vector<isobar> & previous_list);
        // Calculate the integral along the curved secment of pinocchio's head
        complex curved_segment(unsigned int basis_id, double s, std::vector<isobar> & previous_list);

        // Save interpolation of the discontinuity calculated elsewhere into the list of iterations
        inline void save_iteration(basis_grid & grid){ _iterations.push_back(new_iteration( _kinematics, grid, _settings)); };
        inline void skip_iteration(){ _iterations.push_back(new_iteration()); };
        
        // -----------------------------------------------------------------------
        private:

        // Private method only accessible to raw_amplitude
        friend class solver;

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

        // Saved interpolation of phase_shifts
        iterateKT::phase_shift _delta;

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