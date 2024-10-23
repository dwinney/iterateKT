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
#include "utilities.hpp"

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
        {};

        // -----------------------------------------------------------------------
        // Getters for masses

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

        // Special points along the pinnochio path
        inline double A(){ return sth(); };
        inline double B(){ return (M2()-m2())/2; };
        inline double C(){ return pth(); };
        inline double D(){ return rth(); };

        // Analytic continuation of barrier factor along real line
        complex kacser(double s);

        // Bounds of integration in the complex plane
        complex s_plus (double s);
        complex s_minus(double s);

        // Reformulation of the curved area in terms of angular variable phi
        complex phi_plus (double s);
        complex phi_minus(double s);

        // -----------------------------------------------------------------------
        private: 
        
        // Masses
        double _m_parent = 0, _m_daughter = 0;
    };
};

#endif // KINEMATICS_HPP