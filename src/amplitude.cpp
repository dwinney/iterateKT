// Top level class which defines the 1 -> 3 decay process. 
// The template feeds in all the relevant user info for a specific process.
// So far we require all equal mass particles
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "amplitude.hpp"

namespace iterateKT
{
    // -----------------------------------------------------------------------
    // Evaluating an amplitude just sums isobars and their associated prefactors
    // in each channel.

    complex raw_amplitude::evaluate(complex s, complex t, complex u)
    {
        complex result = 0;

        // S_CHANNEL
        for (auto f : _isobars)
        {
            complex term = prefactor_s(f->get_id(), s, t, u);
            if (is_zero(term)) continue;
            result += term * f->evaluate(s);
        };

        // T_CHANNEL
        for (auto f : _isobars)
        {
            complex term = prefactor_t(f->get_id(), s, t, u);
            if (is_zero(term)) continue;
            result += term * f->evaluate(t);
        };

        // U_CHANNEL
        for (auto f : _isobars)
        {
            complex term = prefactor_u(f->get_id(), s, t, u);
            if (is_zero(term)) continue;
            result += term * f->evaluate(u);
        };

        return result;
    };

    // -----------------------------------------------------------------------
    // Calculate partial widths and the integrated width

    // Doubly differential 
    double raw_amplitude::differential_width(double s, double t)
    {
        double u = _kinematics->Sigma() - s - t;

        bool in_physical_region = (std::real(_kinematics->kibble(s, t, u)) >= 0);
        if (!in_physical_region)
        {
            return error("amplitude::differential_width", 
                         "Evaluating outside decay region!", NaN<double>());
        };

        return norm(evaluate(s, t, u))/prefactors()/helicity_factor();
    };

    // Singly differential 
    double raw_amplitude::differential_width(double s)
    {
        using namespace boost::math::quadrature;

        bool in_physical_region = (s >= _kinematics->sth() || s <= _kinematics->pth());
        if (!in_physical_region)
        {
            return error("amplitude::differential_width", 
                         "Evaluating outside decay region!", NaN<double>());
        };

        auto fdx = [&](double t)
        {
            double u = _kinematics->Sigma() - s - t;
            return norm(evaluate(s, t, u))/prefactors()/helicity_factor();
        };

        // Limits are purely real in the decay region
        double min = std::real(_kinematics->t_minus(s));
        double max = std::real(_kinematics->t_plus(s));
        return gauss_kronrod<double,N_GAUSS_ANGULAR>::integrate(fdx, min, max, 0, 1.E-9, NULL);
    };

    // Fully integrated width
    double raw_amplitude::width()
    {
        using namespace boost::math::quadrature;

        auto fdx = [&](double s)
        {
            return differential_width(s);
        };

        double min = _kinematics->sth();
        double max = _kinematics->pth();
        return gauss_kronrod<double,N_GAUSS_ANGULAR>::integrate(fdx, min, max, 0, 1.E-9, NULL);
    };
}; // namespace iterateKT