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

#ifndef KINEMATICS_HPP
#define KINEMATICS_HPP

#include <memory>
#include <Math/Interpolator.h>
#include <TMath.h>
#include "utilities.hpp"
#include "plotter.hpp"
#include "plot2D.hpp" 

namespace iterateKT
{
    // Forward declare amplitude for the typedef below
    class raw_kinematics;

    // Define amplitudes only as pointers
    using kinematics  = std::shared_ptr<raw_kinematics>;

    // This function serves as our "constructor"
    inline kinematics new_kinematics(double m_parent, double m_daughter)
    {
        return std::make_shared<raw_kinematics>(m_parent, m_daughter);
    };

    class raw_kinematics
    {
        // -----------------------------------------------------------------------
        public: 

        raw_kinematics(double m_parent, double m_daughter)
        : _m_parent(m_parent), _m_daughter(m_daughter)
        {
            initialize();
        };

        // -----------------------------------------------------------------------
        // Related to masses

        inline double M(){ return _m_parent;   };
        inline double m(){ return _m_daughter; };
        
        // Masses squared
        inline double m2(){ return m()*m(); };
        inline double M2(){ return M()*M(); };

        
        // Threshold & pseudo-threshold and final state threshold
        inline double sth() { return 4.*m2(); };
        inline double pth() { return (M()-m())*(M()-m()); };
        inline double rth() { return (M()+m())*(M()+m()); };

        // Sum of masses squared
        // s + t + u = Sigma (note no factor of 3!)
        inline double Sigma(){ return M2() + 3*m2(); };
        inline double r(){ return Sigma()/3; };
        inline double s0(){return r(); }; // Equivalent name 

        // Kibble function
        inline complex kibble(complex s, complex t, complex u){ return s*t*u - m2()*pow(M2()-m2(), 2); };

        // Special points along the pinnochio path
        inline double A(){ return sth(); };
        inline double B(){ return (M2()-m2())/2; };
        inline double C(){ return pth(); };
        inline double D(){ return rth(); };

        // -----------------------------------------------------------------------
        // Related to momenta

        // The modulus of 3-momentum for initial or final state in CM frame
        inline complex momentum_initial(complex s){ return csqrt(kallen(s, M2(), m2()))/csqrt(4*s); };
        inline complex momentum_final  (complex s){ return csqrt(kallen(s, m2(), m2()))/csqrt(4*s); };

        // Analytic continuation of barrier factor along real line
        complex kacser   (complex s);
        inline complex kz(complex s, complex t){ return 2*t + s - M2() - 3*m2(); };

        // Kacser with removed singularities at regular thresholds removed
        // xi is the radius of validity to remove the singularities analytically
        complex nu (double s, std::array<double,3> xi = {EPS, EPS, EPS});

        // -----------------------------------------------------------------------
        // Related to boundary of physical regions

        // Kibble function
        inline complex kibble(complex s, complex t, complex u){ return s*t*u - m2()*pow(M2()-m2(), 2); };
        inline bool in_decay_region(double s, double t){ return real(kibble(s,t, Sigma()-s-t)) >= 0; };

        // Bounds of integration in the complex plane
        complex t_plus (double s);
        complex t_minus(double s);
        complex t_curve (double phi);
        complex jacobian(double phi);

        // Reformulation of the curved area in terms of angular variable phi
        double phi_plus (double s);
        double phi_minus(double s);
        double radius(double phi);

        // -----------------------------------------------------------------------
        // Plotting utility for decay region

        // Set up a plot2D in the physical decay region
        plot2D new_dalitz_plot(plotter & pltr, int N = 300);

        // -----------------------------------------------------------------------
        private: 
        
        // Set up an interpolation of the phi contour so that we can take derivatives
        void initialize();
        bool _initialized = false;
        int  _n_interp    = 100;
        ROOT::Math::Interpolator _re_tphi, _im_tphi;
        ROOT::Math::Interpolator _re_jac,  _im_jac;

        // Masses
        double _m_parent = 0, _m_daughter = 0;
    };
};

#endif // KINEMATICS_HPP