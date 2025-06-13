// Amplitudes relevant for the decay of meson with JP = 1-+ into 3pi
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef PI1_AMPLITUDES_HPP
#define PI1_AMPLITUDES_HPP

#include "amplitude.hpp"
#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "GKPY.hpp"

#include"isobars/pi1.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace iterateKT
{ 
    // This defines the full amplitude, i.e. how the isobars are combined
    // Here is where we usually put the isospin combinations etc
    class pi1 : public raw_amplitude
    {
        public: 
        
        // Constructor
        pi1(kinematics kin, std::string id) : raw_amplitude(kin,id)
        {};
        
        // Spin 1 decay so (2j+1) = 3
        inline double combinatorial_factor(){ return 3; };

        // Loop integral
        static inline complex projected_deck(double t, double m3pi2, complex sig)
        {
            using namespace boost::math::quadrature;
            auto f = [&](double x, double y)
            {
                complex a, b, c, d, mu = M_PION;
                a = t;
                b = x*(mu*mu + t - m3pi2) - t;
                c = (1-x)*(1-x)*mu*mu + x*sig;
                d = csqrt(b*b - 4*a*c);
                return (4*a*y - (b+2*a*y)*log(a*y*y+b*y+c) - 2*d*atanh((b+2*a*y)/d))/a;
            };
            auto integrand = [&](double x){ return f(x, 1-x) - f(x, 0); };
            return gauss_kronrod<double,61>::integrate(integrand, 0, 1, 0, 1.E-9, NULL);
        };
        
        // Assuming a pi- pi- pi+ decay and only P-waves
        // s = (pi- + pi+)^2 
        // t = (pi- + pi+)^2
        // u = (pi- + pi-)^2
        inline complex prefactor_s(id iso_id, complex s, complex t, complex u){ return csqrt(_kinematics->kibble(s,t,u)); };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u){ return - prefactor_s(iso_id, t, s, u); };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u){ return 0.; };
    };
}; // namespace iterateKT 

#endif // PI1_AMPLITUDES_HPP