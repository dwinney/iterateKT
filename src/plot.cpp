// Reworking of the jpacStyle as internal library special for jpacPhoto
// This defines a single plot object, which is what actually prints out the files
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "plot.hpp"
#include "colors.hpp"

namespace iterateKT
{
    // ---------------------------------------------------------------------------
    // Outputting function, generates plot and saves it to file
    void plot::save()
    {
        if (_entries.empty() ) 
        {
            warning("plot::save()", "No entries added! Returning...");
            return;
        };

        _canvas->cd();
        
        // Logscale settings are global canvas settings so do that first
        _canvas->SetLogx(_xlog); _canvas->SetLogy(_ylog);
        
        // Call all the graph methods and draw them onto canvas
        draw();

        // Draw the canvas
        _canvas->Draw();

        // and print to file
        _canvas->Print(_filename.c_str());
    };

    void plot::draw()
    {
        bool add_padding = _customyrange && std::abs(_ybounds[0]) < 0.09 && std::abs(_ybounds[1]) < 0.09;
        if (add_padding) 
        {
            auto style = gROOT->GetStyle("jpacStyle");
            style->SetTitleOffset(1.8, "y");
            _canvas->UseCurrentStyle();
        };

        TMultiGraph * mg = new TMultiGraph("mg", "mg");

        // Set up the axes by grabbing them from the first entry
        std::string labels = ";" + _xlabel + ";" + _ylabel;
        mg->SetTitle(labels.c_str());

        // Draw the first entry
        // need to parse if its a curve or data points

        // Set up legend
        _legendxoffset  = 0.3;
        _legendyoffset  = _legendyscale*(_Nlegend + _addheader);
        auto legend = new TLegend(_legendxcord,  _legendycord, 
                                  _legendxcord + _legendxoffset, 
                                  _legendycord + _legendyoffset);
        legend->SetFillStyle(0); // Make legend transparent

        // IF theres a custom header add it
        if (_addheader) legend->SetHeader(("  " + _header).c_str(), "L");

        for (auto entry : _entries)
        {
            mg->Add(entry._graph, entry._style._draw_opt.c_str());

            if (entry._style._add_to_legend)
            {
                legend->AddEntry(entry._graph, entry._style._label.c_str(), entry._style._draw_opt.c_str());
            };
        };
        mg->Draw("same");

        if (_add_logo) add_logo();

        if (_prelim) add_watermark();

        if (_addlegend) legend->Draw();
        if (!_addlegend && _addheader) { legend->Clear(); legend->SetHeader(("  " + _header).c_str(), "L"); legend->Draw();}

        mg->GetXaxis()->CenterTitle(true);
        mg->GetYaxis()->CenterTitle(true);
        mg->GetYaxis()->SetMaxDigits(2);
        _canvas->Modified();

        double ylow, yhigh, xlow, xhigh;
        if (_customxrange)
        {
            mg->GetXaxis()->SetLimits(   _xbounds[0], _xbounds[1]);
            xlow = _xbounds[0]; xhigh = _xbounds[1];
        };

        if (_customyrange)
        {
                        mg->GetYaxis()->SetRangeUser(_ybounds[0], _ybounds[1]);
            mg->SetMinimum(_ybounds[0]);
            mg->SetMaximum(_ybounds[1]);
            ylow = _ybounds[0]; yhigh = _ybounds[1];
        }
        else
        {
            ylow  = mg->GetYaxis()->GetXmin();
            yhigh = mg->GetYaxis()->GetXmax();
        };

        double xmin = mg->GetXaxis()->GetXmin();
        double xmax = mg->GetXaxis()->GetXmax();
        
        for (auto vline : _vlines)
        {
            auto vert = new TLine(vline._value, ylow, vline._value, yhigh);
            vert->SetLineWidth(0.7*_scale*plot_entry::_default_linewidth);
            vert->SetLineStyle(vline._linestyle);
            vert->SetLineColorAlpha(vline._color, 0.7);
            vert->Draw();
        }

        for (auto hline : _hlines)
        {
            auto hori = new TLine(xlow, hline._value, xhigh, hline._value);
            hori->SetLineWidth(0.7*_scale*plot_entry::_default_linewidth);
            hori->SetLineStyle(hline._linestyle);
            hori->SetLineColorAlpha(hline._color, 0.7);
            hori->Draw();
        }

        for (auto shade : _shaded)
        {
            double bmin = (shade._xmin < xmin) ? xmin : shade._xmin;
            double bmax = (shade._xmax > xmax) ? xmax : shade._xmax;

            TBox *b = new TBox(bmin, ylow, bmax, yhigh); 
            b->SetFillColorAlpha(shade._color, 0.1); 
            b->SetFillStyle(shade._style);
            b->Draw();
        }

        // Remove padding at the end
        if (add_padding) 
        {
            auto style = gROOT->GetStyle("jpacStyle");
            style->SetTitleOffset(1.5, "y");
            _canvas->UseCurrentStyle();
        };
    };

    // ---------------------------------------------------------------------------
    // Convert data_set and amplitude easily into plot_entries

    void plot::add_data(std::array<std::vector<double>,2> vx, std::array<std::vector<double>,2> vy, jpacColor color)
    {
        double *x, *y, *dx, *dy;
        x  = &(vx[0][0]);        y  = &(vy[0][0]);
        dx = &(vx[1][0]);       dy  = &(vy[1][0]);
        TGraph *graph = new TGraphErrors(vx[0].size(), x, y, dx, dy);

        entry_style style;
        style._style = 20 + _Ndata;
        style._color = color;
        style._draw_opt = "P";
        style._add_to_legend = false;

        _Ndata++;
        _Nlegend++;

        _entries.push_front(plot_entry(graph, style, true));
    };

    void plot::add_data(std::vector<double> vx, std::array<std::vector<double>,2> vy, jpacColor color)
    {
        double *x, *y, *dy;
        x  = &(vx[0]);        y  = &(vy[0][0]);
                             dy  = &(vy[1][0]);
        TGraph *graph = new TGraphErrors(vx.size(), x, y, nullptr, dy);

        entry_style style;
        style._style = 20 + _Ndata;
        style._color = color;
        style._draw_opt = "P";
        style._add_to_legend = false;

        _Ndata++;
        _Nlegend++;

        _entries.push_front(plot_entry(graph, style, true));
    };

    // -----------------------------------------------------------------------
    // Add different curves from amplitudes / functions

    // Add a curve from raw vectors of (x, fx) pairs
    // this is the most reminicent of jpacStyle
    // color is automatically cycled through the jpacColors
    void plot::add_curve(std::vector<double> x, std::vector<double> fx, entry_style style)
    {
        if (style._add_to_legend) _Nlegend++;
        TGraph *g = new TGraph(x.size(), &(x[0]), &(fx[0]));
        _entries.push_back(plot_entry(g, style, false));
    };

    void plot::add_curve(std::vector<double> x, std::vector<double> fx, std::string id)
    {
        _Ncurve++;
        entry_style style;
        style._color = JPACCOLORS[_Ncurve];
        style._style = kSolid;
        style._label = id;
        style._add_to_legend = (id != "");
        add_curve(x, fx, style);
    };

    // Take in a lambda an evaluation range to get the vectors
    std::array<std::vector<double>,2> plot::add_curve(std::array<double,2> bounds, std::function<double(double)> F, entry_style style)
    {
        double step = (bounds[1] - bounds[0]) / double(_Npoints);

        std::vector<double> x, fx;
        for (int n = 0; n < _Npoints; n++)
        {
            double xs  = (!_xlog) ? bounds[0] + double(n) * (bounds[1] - bounds[0]) / double(_Npoints-1)
                                  : exp(log(bounds[0]) + double(n) * (log(bounds[1]) - log(bounds[0])) / double(_Npoints-1));
            double fxs = F(xs);

            x.push_back(xs);
            fx.push_back(fxs);

            if (_print) print(xs, fxs);
        };

        add_curve(x, fx, style);

        return std::array<std::vector<double>,2>({x, fx});
    };

    std::array<std::vector<double>,2> plot::add_curve(std::array<double,2> bounds, std::function<double(double)> F, std::string id)
    {
        _Ncurve++;
        entry_style style;
        style._color = JPACCOLORS[_Ncurve];
        style._style = kSolid;
        style._label = id;
        style._add_to_legend = (id != "");

        return add_curve(bounds, F, style);
    };
    
    void plot::add_dashed(std::vector<double> x, std::vector<double> fx)
    {
        entry_style style;
        style._color = JPACCOLORS[_Ncurve];
        style._style = kDashed;
        style._add_to_legend = false;

        TGraph *g = new TGraph(x.size(), &(x[0]), &(fx[0]));
        _entries.push_back(plot_entry(g, style, false));
    };

     std::array<std::vector<double>,2> plot::add_dashed(std::array<double,2> bounds, std::function<double(double)> F)
    {
        double step = (bounds[1] - bounds[0]) / double(_Npoints);

        std::vector<double> x, fx;
        for (int n = 0; n < _Npoints; n++)
        {
            double xs  = (!_xlog) ? bounds[0] + double(n) * (bounds[1] - bounds[0]) / double(_Npoints-1)
                                  : exp(log(bounds[0]) + double(n) * (log(bounds[1]) - log(bounds[0])) / double(_Npoints-1));
            double fxs = F(xs);

            x.push_back(xs);
            fx.push_back(fxs);
            if (_print) print(xs, fxs);
        };

        add_dashed(x, fx);
        return std::array<std::vector<double>,2>({x, fx});
    };

    // -----------------------------------------------------------------------
    // Add an error band

    void plot::add_band(std::vector<double> x, std::array<std::vector<double>,2> band, int fill)
    {
        std::vector<double> higher, lower;
        lower = band[0]; higher = band[1];
        
        int   N = x.size();
        auto  y = (higher + lower) / 2;
        auto ey = (higher - lower) / 2;

        TGraph *graph = new TGraphErrors(N, &(x[0]), &(y[0]), NULL, &(ey[0]));

        jpacColor color = JPACCOLORS[_Ncurve];
        graph->SetFillColorAlpha(+color, 0.25);
        graph->SetFillStyle(fill);

        entry_style style;
        style._color = color;
        style._draw_opt = "3";
        style._add_to_legend = false;
        _entries.push_front(plot_entry(graph, style, false));
    };
};