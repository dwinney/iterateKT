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

    // -----------------------------------------------------------------------
    // Automate making plots of the amplitude

    // Instead of outputting an array, we output a vector to make it 
    // compatible right away with plotter.combine
    std::vector<plot2D> raw_amplitude::make_plots(plotter & pltr, int N)
    {
        double smin  = _kinematics->sth();
        double smax  = _kinematics->pth();
        double sigma = _kinematics->Sigma();

        std::vector<double> s, t, re, im, abs; 
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
                abs.push_back( norm(ampij) );
            };
        };

        plot2D p_re = _kinematics->new_dalitz_plot(pltr);
        p_re.set_data({s,t,re});
        p_re.set_title("Re#kern[0.2]{(}#it{A})");
        p_re.set_labels("#it{s}", "#it{t}");


        plot2D p_im = _kinematics->new_dalitz_plot(pltr);
        p_im.set_data({s,t,im});
        p_im.set_title("Im#kern[0.2]{(}#it{A})");
        p_im.set_labels("#it{s}", "#it{t}");

        plot2D p_abs = _kinematics->new_dalitz_plot(pltr);
        p_abs.set_data({s,t,abs});
        p_abs.set_title("|#kern[0.3]{#it{A}}|#kern[0.2]{^{2}}");
        p_abs.set_labels("#it{s}", "#it{t}");

        return {p_re, p_im, p_abs};
    };

}; // namespace iterateKT