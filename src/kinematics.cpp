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
    complex raw_kinematics::kacser(complex s)
    {
        complex kappa = csqrt(1-sth()/s)*csqrt(pth()-s)*csqrt(rth()-s);
        if (std::real(s) < sth() || !is_zero(std::imag(s))) return kappa;
        
        // Handle the continuation piecewise
        int region = (std::real(s) >= pth()) + (std::real(s) > rth());
        switch (region)
        {
            case 0: return +  std::abs(kappa);
            case 1: return +I*std::abs(kappa);
            case 2: return -  std::abs(kappa);
        };
        return NaN<complex>();
    };

    // Bounds of integration
    complex raw_kinematics::t_plus(double s)
    {
        return (Sigma() - s + kacser(s))/2;
    };
    complex raw_kinematics::t_minus(double s)
    {
        return (Sigma() - s - kacser(s))/2;
    };

    // -----------------------------------------------------------------------
    // These methods are related to the curved section of pinocchios head

    complex raw_kinematics::t_curve(double phi)
    {
        if (!_initialized) return radius(phi)*exp(I*phi);
        return _re_tphi.Eval(phi)+I*_im_tphi.Eval(phi);
    };

    complex raw_kinematics::jacobian(double phi)
    {
        if (!_initialized) initialize();
        return _re_tphi.Deriv(phi)+I*_im_tphi.Deriv(phi);
    };

    void raw_kinematics::initialize()
    {
        std::vector<double> phi, re, im;
        for (int i = 0; i < _n_interp; i++)
        {
            double phi_i = (2*PI)*double(i)/double(_n_interp-1);
            complex tprime = t_curve(phi_i);

            phi.push_back(phi_i);
            re.push_back( std::real(tprime) );
            im.push_back( std::imag(tprime) );
        };

        _re_tphi.SetData(phi, re);
        _im_tphi.SetData(phi, im);
        _initialized = true;
    };

    // Bounds of integration in terms of phi
    double raw_kinematics::phi_plus(double s)
    {
        if (s >= rth() || s <= pth()) return error("kinematics::phi_plus", 
                                                   "Outside egg region!", NaN<double>());
        double cosine = (Sigma()-s)*sqrt(s)/(M2()-m2())/m();
        return acos(cosine);
    };
    double raw_kinematics::phi_minus(double s)
    {
        return 2*PI - phi_plus(s);
    };

    // The radius of integratoin path t along the curved section
    double raw_kinematics::theta(double phi)
    {
        double arg = cos(phi)*m()*(M2()-m2());
        arg *= arg;
        arg *= 2/pow(Sigma()/3, 3.);
        return acos(arg-1)/3;
    };

    double raw_kinematics::radius(double phi)
    {
        double cosphi   = cos(phi);
        if (is_zero(cosphi)) return (M2() - m2())*m()/Sigma();
        if (cosphi < 0) return (Sigma()/3)*(1-2*cos(theta(phi)))/2/cosphi;
        else return (Sigma()/3)*(1+2*cos(theta(phi)+PI/3.))/2/cosphi;
    };  

}; // namespace iterateKT