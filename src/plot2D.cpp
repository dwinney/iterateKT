// Reworking of the jpacStyle as internal library special for jpacPhoto
// This defines a single plot object, which is what actually prints out the files
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "plot2D.hpp"

namespace iterateKT
{
    void plot2D::save()
    {
        _canvas->cd();

        // Call all the graph methods and draw them onto canvas
        draw();

        // Draw the canvas
        _canvas->Draw();

        // and print to file
        _canvas->Print(_filename.c_str());
    };

    void plot2D::draw()
    {
        TGraph2D* graph = new TGraph2D(_data[0].size(), &(_data[0][0]), &(_data[1][0]), &(_data[2][0]));
        graph->SetName(_filename.c_str());

        graph->SetTitle(_title.c_str());

        TAxis* xAxis = graph->GetHistogram()->GetXaxis();
        xAxis->SetTitle(_xlabel.c_str());
        xAxis->CenterTitle(true);

        TAxis* yAxis = graph->GetHistogram()->GetYaxis();
        yAxis->SetTitle(_ylabel.c_str());
        yAxis->CenterTitle(true);
     
        if (_custom_ranges)
        {
            xAxis->SetLimits(_xbounds[0], _xbounds[1]);
            yAxis->SetLimits(_ybounds[0], _ybounds[1]);
        };

        if (_custom_region)
        {
            // Make an exclusion of the dalitz region
            int Nreg = _region[0].size();
            TCutG *cut = new TCutG("cut", Nreg);
            cut->SetVarX("y");
            cut->SetVarY("x");

            TGraph *outline = new TGraph(Nreg);
            outline->SetLineWidth(2);
            for (int i = 0; i < Nreg; i++) 
            {
                cut->SetPoint(i, _region[0][i], _region[1][i]);
                outline->SetPoint(i, _region[0][i], _region[1][i]);
            };
            graph->Draw("COLZ0 [cut]");
            outline->Draw("L SAME");      
            return;
        };

        graph->Draw("COLZ0");
    };
};