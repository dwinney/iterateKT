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
    // Momentum barrier factors

    // Full kacser function = 4p*q
    complex raw_kinematics::kacser(complex s)
    {
        complex kappa = csqrt(1-sth()/s)*csqrt(pth()-s)*csqrt(rth()-s);
        if (real(s) < sth() || !is_zero(imag(s))) return kappa;
        
        // Handle the continuation piecewise
        int region = (real(s) >= pth()) + (real(s) > rth());
        switch (region)
        {
            case 0: return +  abs(kappa);
            case 1: return +I*abs(kappa);
            case 2: return -  abs(kappa);
        };
        return NaN<complex>();
    };

    // Kacser with removed singularities at regular thresholds removed
    // xi is the radius of validity to remove the singularities analytically
    complex raw_kinematics::nu(double s, std::array<double,3> xi)
    {
        double Sth = sth();
        double Rth = rth();

        int region; 
        region = (s >= Sth + xi[0])  // 1
               + (s >= Rth - xi[2])  // 2
               + (s >= Rth)          // 3
               + (s >= Rth + xi[2]); // 4

        switch (region)
        {
            case 0: return sqrt(Rth - s)/sqrt(s);
            case 1: return sqrt(Rth - s)*sqrt(1 - Sth/s);
            case 2: return   sqrt(1 - Sth/s);
            case 3: return I*sqrt(1 - Sth/s);
            case 4: return I*sqrt(s - Rth)*sqrt(1 - Sth/s);
            default : return NaN<double>();
        };
    };

    // -----------------------------------------------------------------------
    // Bounds of angular integration

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
        return _re_jac.Eval(phi)+I*_im_jac.Eval(phi);
    };

    void raw_kinematics::initialize()
    {
        std::vector<double> phi, re, im;
        for (int i = 0; i < _n_interp; i++)
        {
            double phi_i = (2*PI)*double(i)/double(_n_interp-1);
            complex tprime = t_curve(phi_i);

            phi.push_back(phi_i);
            re.push_back( real(tprime) );
            im.push_back( imag(tprime) );
        };

        _re_tphi.SetData(phi, re);
        _im_tphi.SetData(phi, im);

        phi.clear(); re.clear(); im.clear();
        for (int i = 0; i < _n_interp; i++)
        {
            double phi_i = (2*PI)*double(i)/double(_n_interp-1);

            phi.push_back(phi_i);
            re.push_back( _re_tphi.Deriv(phi_i) );
            im.push_back( _im_tphi.Deriv(phi_i) );
        };

        _re_jac.SetData(phi, re);
        _im_jac.SetData(phi, im);
        _initialized = true;
    };

    // Bounds of integration in terms of phi
    double raw_kinematics::phi_plus(double s)
    {
        if (s > rth() || s < pth()) return error("kinematics::phi_plus", 
                                                   "Outside egg region!", NaN<double>());
        double cosine = (Sigma()-s)*sqrt(s)/(M2()-m2())/m()/2;
        if (!are_equal(abs(cosine), 1)) return TMath::ACos(cosine);
        return (cosine > 0) ? 0 : PI;
    };
    double raw_kinematics::phi_minus(double s)
    {
        return 2*PI - phi_plus(s);
    };

    double raw_kinematics::radius(double phi)
    {
        double cosphi   = cos(phi);
        
        // Calculate angle in theta
        double arg = cos(phi)*m()*(M2()-m2());
        arg *= arg;
        arg *= 2/pow(Sigma()/3, 3.);
        double theta = acos(arg-1)/3;

        if (is_zero(cosphi)) return (M2() - m2())*m()/Sigma();
        if (cosphi < 0) return (Sigma()/3)*(1-2*cos(theta))/2/cosphi;
        else return (Sigma()/3)*(1+2*cos(theta+PI/3.))/2/cosphi;
    };  

    // -----------------------------------------------------------------------
    plot2D raw_kinematics::new_dalitz_plot(plotter & pltr, int N)
    {
        std::array<std::vector<double>,2> dalitz_boundary;
        // Go clockwise around the dalitz region, first tmax
        for (int i = 0; i <= N; i++)
        {
            double s = sth() + (pth() - sth())*i/double(N);
            double t = real(t_plus(s));
    
            dalitz_boundary[0].push_back(s); dalitz_boundary[1].push_back(t);
        };
        // then tmin going backwards in s
        for (int i = 0; i <= N; i++)
        {
            double s = pth() + (sth() - pth())*i/double(N);
            double t = real(t_minus(s));
    
            dalitz_boundary[0].push_back(s); dalitz_boundary[1].push_back(t);
        };

        plot2D p = pltr.new_plot2D();
        p.set_region(dalitz_boundary);
        return p;
    };

}; // namespace iterateKT