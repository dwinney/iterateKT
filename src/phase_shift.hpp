// Read in a phase shift from a file, interpolate it, and match it to some asymptotic
// 
// ---------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef PHASE_SHIFT_HPP
#define PHASE_SHIFT_HPP

#include "constants.hpp"
#include "data_set.hpp"
#include <Math/Interpolator.h>
#include <tuple>

namespace iterateKT
{
    // Short-cut for passing arguments in a single structure
    using phase_args = std::tuple<std::string,double,int,double,double>; 

    class phase_shift
    {
        public:

        phase_shift(): _error(true) {};

        // constructor takes in the file name, matching energy, and integrer of pi 
        phase_shift(std::string file, double lam2, uint k, double tau, double a) 
        : _error(false), _match(lam2), _k(k), _tau(tau), _a(a)
        { interpolate(file); };

        phase_shift(phase_args info): _error(false),         _match(std::get<1>(info)), 
                                      _k(std::get<2>(info)), _tau  (std::get<3>(info)), _a(std::get<4>(info))
        { interpolate(std::get<0>(info)); };
        
        inline double operator()(double s)
        {
            if (_error)      return error("phase_shift", "No interpolation specified!", 0.);
            if (s <= _sth)   return 0.;
            if (s <= _match) return _delta.Eval(s);
            return asymptotic(s);
        };

        inline void set_info(phase_args info)
        {
            _error = false; 
            interpolate(std::get<0>(info));
            _match = std::get<1>(info);
            _k     = std::get<2>(info);
            _tau   = std::get<3>(info);
            _a     = std::get<4>(info);
        };

        private:

        bool   _error = true;
        uint   _k;           // Multiple of pi to extrapolate at infinity
        double _match, _sth; // Cutoff and threshold
        double _tau, _a;     // Decay rate and power in the asymptotic matching
        ROOT::Math::Interpolator _delta; 

        inline void interpolate(std::string file)
        {   
            // Assume data is in two columns
            auto data = import_data<2>(file);
            check<2>(data, file);
            _sth = data[0][0]; 
            _delta.SetData(data[0], data[1]);
        };

        inline double asymptotic(double s)
        {
            // This can be any function so long as it 
            // and its first derivative vanish at s = _match;
            double b = pow(_tau*(s-_match),_a);
            return _delta.Eval(_match)*exp(-b)+(1-exp(-b))*_k*PI;
        };
    };
};
#endif