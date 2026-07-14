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

    complex raw_amplitude::evaluate_precomputed(double s, double t, double u)
    {
        complex result = 0;

        // S_CHANNEL
        for (auto f : _isobars)
        {
            complex term = prefactor_s(f->get_id(), s+IEPS, t+IEPS, u+IEPS);
            if (is_zero(term)) continue;
            result += term * f->evaluate_precomputed(s);
        };

        // T_CHANNEL
        for (auto f : _isobars)
        {
            complex term = prefactor_t(f->get_id(), s+IEPS, t+IEPS, u+IEPS);
            if (is_zero(term)) continue;
            result += term * f->evaluate_precomputed(t);
        };

        // U_CHANNEL
        for (auto f : _isobars)
        {
            complex term = prefactor_u(f->get_id(), s+IEPS, t+IEPS, u+IEPS);
            if (is_zero(term)) continue;
            result += term * f->evaluate_precomputed(u);
        };

        return result;
    };

    // -----------------------------------------------------------------------
    // Calculate partial widths and the integrated width

    // Doubly differential 
    double raw_amplitude::differential_width(double s, double t)
    {
        bool in_physical_region = get_kinematics()->in_decay_region(s,t);
        if (!in_physical_region)
        {
            return error("amplitude::differential_width", 
                         "Evaluating outside decay region!", NaN<double>());
        };

        return norm(evaluate_in_dalitz(s, t))/prefactors();
    };

    // Singly differential 
    double raw_amplitude::differential_width(double s)
    {
        using namespace boost::math::quadrature;

        bool in_physical_region = (s >= get_kinematics()->sth() || s <= get_kinematics()->pth());
        if (!in_physical_region)
        {
            return error("amplitude::differential_width", 
                         "Evaluating outside decay region!", NaN<double>());
        };

        auto fdx = [&](double t)
        {
            return differential_width(s, t);
        };

        // Limits are purely real in the decay region
        double min = real(get_kinematics()->t_minus(s));
        double max = real(get_kinematics()->t_plus(s));
        return gauss_kronrod<double,61>::integrate(fdx, min, max, 0, 1.E-9, NULL);
    };

    // Fully integrated width
    double raw_amplitude::width()
    {
        using namespace boost::math::quadrature;

        auto fdx = [&](double s){ return differential_width(s); };

        double min = get_kinematics()->sth();
        double max = get_kinematics()->pth();
        return gauss_kronrod<double,61>::integrate(fdx, min, max, 0, 1.E-9, NULL);
    };

    // -----------------------------------------------------------------------
    // Automate making plots of the amplitude

    // Instead of outputting an array, we output a vector to make it 
    // compatible right away with plotter.combine
    std::vector<plot2D> raw_amplitude::plot_ReIm(plotter & pltr, std::string units, int N)
    {
        double smin  = get_kinematics()->sth();
        double smax  = get_kinematics()->pth();
        double sigma = get_kinematics()->Sigma();

        std::vector<double> s, t, re, im; 
        for (int i = 0; i < N; i++)
        {
            double si = smin+(smax-smin)*i/double(N-1);

            double tmin = real(get_kinematics()->t_minus(si));
            double tmax = real(get_kinematics()->t_plus (si));
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

        plot2D p_re = get_kinematics()->new_dalitz_plot(pltr);
        p_re.set_data({s,t,re});
        p_re.set_title("Re#kern[0.2]{(}#it{A})");
        p_re.set_labels(xlabel, ylabel);


        plot2D p_im = get_kinematics()->new_dalitz_plot(pltr);
        p_im.set_data({s,t,im});
        p_im.set_title("Im#kern[0.2]{(}#it{A})");
        p_im.set_labels(xlabel, ylabel);

        return {p_re, p_im};
    };

    // Plot |A|
    plot2D raw_amplitude::plot_dalitz(plotter & pltr, std::string units, int N)
    {
        double smin  = get_kinematics()->sth();
        double smax  = get_kinematics()->pth();
        double sigma = get_kinematics()->Sigma();

        std::vector<double> s, t, absA; 
        for (int i = 0; i < N; i++)
        {
            double si = smin+(smax-smin)*i/double(N-1);

            double tmin = real(get_kinematics()->t_minus(si));
            double tmax = real(get_kinematics()->t_plus (si));
            for (int j = 0; j < N; j++)
            {
                double tij = tmin+(tmax-tmin)*j/double(N-1);
                
                complex ampij = evaluate_precomputed(si, tij, sigma - si - tij);

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

        plot2D p = get_kinematics()->new_dalitz_plot(pltr);
        p.set_data({s,t,absA});
        p.set_labels(xlabel, ylabel);

        return p;
    };

    // -----------------------------------------------------------------------
    // Calculate Daltiz plot parameters from amplitude
    // By default we use the conventions and notation of the PDG 
    // see ``Dalitz Plot Parameters for K -> 3pi decays" in RPP
    // Output in order {g, h, j, k, f}

    // The optional arguments are for the center of the dalitz plot and normalization:
    // X = (t - u)/m[0]
    // Y = (s - r)/m[1]
    
    std::array<double,5> raw_amplitude::get_dalitz_parameters(double e, double s0, std::array<double,2> m)
    {
        double N  = norm(evaluate(s0,s0,s0));

        // Rename our function for readibility
        auto F  = [this,N,s0](double s, double t){ return norm(evaluate(s,t,get_kinematics()->Sigma()-s-t))/N; };
        auto Fs = [this,F,s0](double s){ return F(s,s0); };
        auto Ft = [this,F,s0](double t){ return F(s0,t); };

        double dFds, dFdt, d2Fd2s, d2Fd2t, d2Fdsdt;
        
        // We just use a (4-point) central finite difference 
        // since these are assumed to be well behaved and fairly smooth
        // derivatives of our F in terms of s and t
        dFds    = central_difference_derivative<double>(1, Fs, s0, e);
        dFdt    = central_difference_derivative<double>(1, Ft, s0, e);

        // 2nd Derivatives
        d2Fd2s  = central_difference_derivative<double>(2, Fs, s0, e);
        d2Fd2t  = central_difference_derivative<double>(2, Ft, s0, e);
        d2Fdsdt = mixed_central_derivatives<double>(F, {s0, s0}, e);

        // derivatives of s and t with respect to X and Y
        // dsdX = 0,  dtdX = cx
        // -2*dsdY =  dtdY = cy
        double cx = +m[0]/2, cy = -m[1]/2; 

        // derivatives of F in terms of X and Y
        double dFdX, dFdY, d2Fd2X, d2Fd2Y, d2FdXdY;
        dFdX    = cx* dFdt;
        dFdY    = cy*(dFdt-2*dFds);
        d2Fd2X  = cx*cx* d2Fd2t;
        d2FdXdY = cx*cy*(d2Fd2t - 2*d2Fdsdt);
        d2Fd2Y  = cy*cy*(d2Fd2t - 4*d2Fdsdt + 4*d2Fd2s); 
        
        // Assmble our outputs
        double g, h, j, f, k;
        g = dFdY;
        h = d2Fd2Y/2;
        j = dFdX;
        k = d2Fd2X/2;
        f = d2FdXdY;

        return {g, h, j, k, f};
    };
}; // namespace iterateKT