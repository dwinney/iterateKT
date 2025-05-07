// Simple container classes to pass around interpolations of basis amplitudes 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef BASIS_HPP
#define BASIS_HPP

#include <functional>
#include "utilities.hpp"

namespace iterateKT
{
    // Forward declare the id type
    enum class id : unsigned int;

    class raw_subtractions;
    using subtractions = std::shared_ptr<raw_subtractions>;

    // When handling subtractions we need several pieces of information which can be passed around 
    // with this object
    class raw_subtractions
    {
        public: 

        raw_subtractions(){};
        
        inline unsigned int N_basis()                   { return _ids.size(); };
        inline id       get_id(uint i)                  { return _ids[i]; };
        inline complex  driving_term(uint i, complex x) { return _drivers[i](x); };
        inline complex  get_par(uint i)                 
        {
            if (_values.size() != N_basis())
            {
                if (!_error_thrown)
                {
                    _error_thrown = true;
                    return error("Subtraction coefficients not set!", NaN<complex>());
                };
                return NaN<complex>();                 
            };
            _error_thrown = false;
            return _values[i]; 
        };

        private: 

        friend class solver;
        friend class raw_amplitude;
        
        // TODO: Probably better way to organize this maybe tuples idk
        bool _error_thrown = false;
        std::vector<id>           _ids;                         // The id of the isobar this subtraction coeff appears in
        std::vector<std::function<complex(complex)>> _drivers;  // Driving terms. in most cases these are the polynomials of fixed integer order
        std::vector<std::string>  _names;                       // Optional to give each parameter a name
        std::vector<complex>      _values;                      // Its actual value determined by fit or matching
    };

    struct basis_grid
    {
        basis_grid(){};

        inline unsigned int N_basis(){ return _re_list.size(); };
        
        int _n_singularity;
        std::vector<double> _s_list, _s_around_pth;
        std::vector<std::vector<double>> _re_list, _im_list;
    };
};

#endif //BASIS_HPP