// General class for all things kinematics.
// Basically anything that only depends on the masses involved.
// This will allow us to pass along all kinematics to the full amplitude as well
// as all isobars individually without much copy/pasting
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "kinematics.hpp"

namespace iterateKT
{
    // -----------------------------------------------------------------------
    complex raw_kinematics::kacser(double s)
    {
        if (s < sth()) return error("kinematics::kacser", "Trying to evaluate Kacser below threshold!", NaN<complex>());
        double kappa = std::abs(csqrt(1-sth()/s)*csqrt(pth()-s)*csqrt(rth()-s));
        
        // Handle the continuation piecewise
        int region = (s >= pth()) + (s > rth());
        switch (region)
        {
            case 0: return +  kappa;
            case 1: return +I*kappa;
            case 2: return -  kappa;
        };
    };

    // Bounds of integration
    complex raw_kinematics::s_plus(double s)
    {
        return (Sigma() - s + kacser(s))/2;
    };
    complex raw_kinematics::s_minus(double s)
    {
        return (Sigma() - s - kacser(s))/2;
    };

    // Bounds of integration in terms of phi
    complex raw_kinematics::phi_plus(double s)
    {
        double cosine = (Sigma()-s)*sqrt(s)/(M2()-m2())/m();
        return acos(cosine);
    };
    complex raw_kinematics::phi_minus(double s)
    {
        return 2*PI - phi_plus(s);
    };

}; // namespace iterateKT