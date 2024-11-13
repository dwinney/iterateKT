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

namespace iterateKT
{
    class phase_shift
    {
        public:

        // constructor takes in the file name, matching energy, and integrer of pi 
        phase_shift(std::string file, double lam2, uint k) : _match(lam2), _k(k)
        {
            interpolate(file, lam2);
        };
        
        inline double operator()(double s)
        {
            if (s <= _sth)   return 0.;
            if (s <= _match) return _delta.Eval(s);
            return asymptotic(s);
        };

        private:

        uint    _k;    // Multiple of pi to extrapolate at infinity
        double _match; // Cutoff
        double _sth, _a = 0., _b = 0.;
        ROOT::Math::Interpolator _delta; 

        inline void interpolate(std::string file, double lam2)
        {   
            // Assume data is in two columns
            auto data = import_data<2>("/physics/phase_shifts/"+file);
            check<2>(data, file);

            _sth = data[0][0]; 
            if (data[0].back() < lam2 || lam2 < _sth) warning("phase_shift", "Cutoff outside interpolation range!");

            // Interpolate
            _delta.SetData(data[0], data[1]);

            // Calculate parameters for the matching
            double d = _delta.Eval(_match), dp = _delta.Deriv(_match);

            if (_k == 0) return;

            _a = pow(_k*PI-d, 2)/_match/dp;
            _b =    (_k*PI-d)   /_match/dp - 1;
            if (_b < -1) warning("phase_shift", "Asymptotic matching exhibits pole!");
        };

        inline double asymptotic(double s)
        {
            return _k*PI - _a/(_b+s/_match);
        };
    };
};
#endif