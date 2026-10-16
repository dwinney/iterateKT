// Top (abstract) level class which defines an X -> π transition form factor.
// We take an X -> 3π amplitude object and specify an isobar.
// This defines the type of current being calculated (i.e. P-wave isobar -> vector FF)
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Instituto de Ciencias Nucleares (ICN)
//               Universidad Autónoma Nacional de México (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#ifndef FORM_FACTOR_HPP
#define FORM_FACTOR_HPP

#include <memory>
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "settings.hpp"
#include "isobar.hpp"
#include "utilities.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace iterateKT
{
    // Forward declare for the typedef below
    class raw_form_factor;
    
    // Define form_factor objects only as pointers
    using form_factor = std::shared_ptr<raw_form_factor>;

    // "Constructors"
    template<class A =raw_form_factor>
    inline form_factor new_form_factor(uint nsub, amplitude amp, id proj, settings settings = default_settings())
    {
        auto x = std::make_shared<A>(nsub, amp, proj, settings);
        return std::static_pointer_cast<raw_form_factor>(x);
    };

    class raw_form_factor
    {
        // -----------------------------------------------------------------------
        public:
        
        // A typical form factor will require:
        // - number of subtractions to consider
        // - a decay amplitude (which specifies the Xπ -> ππ amplitude)
        // - the projection id (which specifies the projection i.e. FF with which spin-J)
        raw_form_factor(uint nsub, amplitude amp, id projection, settings settings):
        _n_subtractions(nsub),
        _decay_amplitude(amp), 
        _kinematics(amp->get_kinematics()),
        _settings(settings)
        {
            _projections.push_back(projection);
            _direct_isobars.push_back(amp->get_isobar(projection));
        };

        // OR if multiple isobars contribute to the "direct channel"
        raw_form_factor(uint nsub, amplitude amp, std::vector<id> projections, settings settings):
        _n_subtractions(nsub),
        _decay_amplitude(amp), 
        _projections(projections),
        _kinematics(amp->get_kinematics()),
        _settings(settings)
        {
            for (auto x : projections)
            {
                _direct_isobars.push_back(amp->get_isobar(x));
            };
        };
        
        // -----------------------------------------------------------------------
        // These functions need to be implemented by a user-defined derived class 
        
        // In addition to the decay amplitude, we need to specify the
        // ππ -> J form factor for arbitrary s on the real line
        virtual complex external_current(double s) = 0; 

        // Depending on how the amplitude is normalizes, we may need to
        // add extra factors to match the relevant partial wave projection
        virtual double kinematic_factors(double s) = 0;

        // -----------------------------------------------------------------------
        // With the above specified, the evaluation simply comes down 
        // to evaluating a dispersion relation

        complex evaluate(complex s); 
        
        // Discontinuity across RHC, could be singular
        complex discontinuity(double s);

        // -----------------------------------------------------------------------
        // Parameter settings
        
        void set_parameters(std::vector<complex> pars)
        {
            if (pars.size() != _n_subtractions)
            {
                warning("form_factor", "Wrong number of parameters recieved!");
                return;
            };
            _subtractions = pars;
        };


        // -----------------------------------------------------------------------
        protected:
        
        // Holds all parameters related to integration and expansions and etc.
        settings _settings; 
        
        // Related to subtraction polynomials
        uint                 _n_subtractions; // number of subtractions
        std::vector<complex> _subtractions;   // subtraction coefficients
        
        // Get the kinematics from the decay amplitude
        kinematics _kinematics; 
        
        // Decay amplitude which supplies the Xπ -> ππ partial waves
        amplitude _decay_amplitude;
        
        // This id specifies the "direct channel" isobar and the relevant partial-wave projection
        std::vector<id> _projections; 
        
        // Save a pointer to the direct channel isobar for ease
        std::vector<isobar> _direct_isobars;
        
        // Functions for evaluation
        
        // Evaluation of the non-singular part of the dispersion integral
        complex regular_piece(complex s);
    };
};

#endif