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

        bool in_physical_region = (real(_kinematics->kibble(s, t, u)) >= 0);
        if (!in_physical_region)
        {
            return error("amplitude::differential_width", 
                         "Evaluating outside decay region!", NaN<double>());
        };

        return norm(evaluate(s, t, u))/prefactors();
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
            return norm(evaluate(s, t, u))/prefactors();
        };

        // Limits are purely real in the decay region
        double min = real(_kinematics->t_minus(s));
        double max = real(_kinematics->t_plus(s));
        return gauss_kronrod<double,15>::integrate(fdx, min, max, 0, 1.E-9, NULL);
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
        return gauss_kronrod<double,15>::integrate(fdx, min, max, 0, 1.E-9, NULL);
    };

    // -----------------------------------------------------------------------
    // Automate making plots of the amplitude

    // Instead of outputting an array, we output a vector to make it 
    // compatible right away with plotter.combine
    std::vector<plot2D> raw_amplitude::plot_ReIm(plotter & pltr, std::string units, int N)
    {
        double smin  = _kinematics->sth();
        double smax  = _kinematics->pth();
        double sigma = _kinematics->Sigma();

        std::vector<double> s, t, re, im; 
        for (int i = 0; i < N; i++)
        {
            double si = smin+(smax-smin)*i/double(N-1);

            double tmin = real(_kinematics->t_minus(si));
            double tmax = real(_kinematics->t_plus (si));
            for (int j = 0; j < N; j++)
            {
                double tij = tmin+(tmax-tmin)*j/double(N-1);
                
                complex ampij = evaluate(si, tij, sigma - si - tij);

                s.push_back(si); t.push_back(tij);
                re.push_back(  real(ampij) );
                im.push_back(  imag(ampij) );
            };
        };

        std::string xlabel = "#sigma_{1}";
        std::string ylabel = "#sigma_{2}";

        if (units != "")
        {
            xlabel += " " + units;
            ylabel += " " + units;
        }

        plot2D p_re = _kinematics->new_dalitz_plot(pltr);
        p_re.set_data({s,t,re});
        p_re.set_title("Re#kern[0.2]{(}#it{A})");
        p_re.set_labels(xlabel, ylabel);


        plot2D p_im = _kinematics->new_dalitz_plot(pltr);
        p_im.set_data({s,t,im});
        p_im.set_title("Im#kern[0.2]{(}#it{A})");
        p_im.set_labels(xlabel, ylabel);

        return {p_re, p_im};
    };

    // Plot |A|
    plot2D raw_amplitude::plot_dalitz(plotter & pltr, std::string units, int N)
    {
        double smin  = _kinematics->sth();
        double smax  = _kinematics->pth();
        double sigma = _kinematics->Sigma();

        std::vector<double> s, t, absA; 
        for (int i = 0; i < N; i++)
        {
            double si = smin+(smax-smin)*i/double(N-1);

            double tmin = real(_kinematics->t_minus(si));
            double tmax = real(_kinematics->t_plus (si));
            for (int j = 0; j < N; j++)
            {
                double tij = tmin+(tmax-tmin)*j/double(N-1);
                
                complex ampij = evaluate(si, tij);

                s.push_back(si); t.push_back(tij);
                absA.push_back(  abs(ampij) );
            };
        };

        std::string xlabel = "#sigma_{1}";
        std::string ylabel = "#sigma_{2}";

        if (units != "")
        {
            xlabel += " " + units;
            ylabel += " " + units;
        }

        plot2D p = _kinematics->new_dalitz_plot(pltr);
        p.set_data({s,t,absA});
        p.set_labels(xlabel, ylabel);

        return p;
    };

    // -----------------------------------------------------------------------
    // Calculate Daltiz plot parameters from amplitude
    // We use the conventions and notation of the PDG 
    // see ``Dalitz Plot Parameters for K -> 3pi decays" in RPP
    // Output in order {g, h, j, k, f}
    
    std::array<double,5> raw_amplitude::get_dalitz_parameters(double e, double m2)
    {
        double s0 = _kinematics->s0(); 
        double N  = norm(evaluate(s0,s0));

        // Rename our function for readibility
        auto F  = [this,N,s0](double s, double u){ return norm(evaluate(s,3*s0-s-u))/N; };
        auto Fs = [this,F,s0](double s){ return F(s,s0); };
        auto Fu = [this,F,s0](double u){ return F(s0,u); };

        
        double dFds, dFdu, d2Fd2s, d2Fd2u, d2Fdsdu;
        
        // We just use a (4-point) central finite difference 
        // since these are assumed to be well behaved and fairly smooth
        // derivatives of our F in terms of s and u
        dFds    = central_difference_derivative<double>(1, Fs, s0, e);
        dFdu    = central_difference_derivative<double>(1, Fu, s0, e);

        // 2nd Derivatives
        d2Fd2s  = central_difference_derivative<double>(2, Fs, s0, e);
        d2Fd2u  = central_difference_derivative<double>(2, Fu, s0, e);
        d2Fdsdu = mixed_partial_derivatives<double>(F, {s0, s0}, e);

        // derivatives of s and u with respect to X and Y
        // X = (t - s) /m2
        // Y = (u - x)/m2
        // dsdX = dsdY = -dudY/2 = c 
        double c = -m2/2; 

        // derivatives of F in terms of X and Y
        double dFdX, dFdY, d2Fd2X, d2Fd2Y, d2FdXdY;
        dFdX    = c* dFds;
        dFdY    = c*(dFds-2*dFdu);
        d2Fd2X  = c*c* d2Fd2s;
        d2FdXdY = c*c*(d2Fd2s - 2*d2Fdsdu);
        d2Fd2Y  = c*c*(d2Fd2s - 4*d2Fdsdu + 4*d2Fd2u); 
        
        // Assmble our outputs
        double g, h, j, f, k;
        g = dFdY;
        h = d2Fd2Y/2;
        j = dFdX;
        k = d2Fd2X/2;
        f = d2FdXdY;

        return {g, h, j, k, f};
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