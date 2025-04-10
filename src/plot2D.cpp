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
        // _canvas->Draw();
        _canvas->Update();
        _canvas->Modified();

        // and print to file
        if (_inverted) TColor::InvertPalette();
        _canvas->Print(_filename.c_str());
        if (_inverted) TColor::InvertPalette();
    };

    void plot2D::draw()
    {
        TGraph2D* graph = new TGraph2D(_data[0].size(), &(_data[0][0]), &(_data[1][0]), &(_data[2][0]));
        graph->SetNpx(_nbins); graph->SetNpy(_nbins);
        gStyle->SetNumberContours(_ncontours);
        graph->SetName(_filename.c_str());

        std::string frame_labels = _title + ";" + _xlabel + ";" + _ylabel + ";";
        graph->SetTitle(frame_labels.c_str());
        graph->GetHistogram()->GetXaxis()->CenterTitle(true);
        graph->GetHistogram()->GetYaxis()->CenterTitle(true);
        std::string draw_options = "COLZ0";

        std::string    command  = "gStyle->SetPalette(" + to_string(_palette)+");";
        TExec *ex = new TExec("ex", command.c_str());

        if (_custom_ranges)
        { 
            auto frame = gPad->DrawFrame(_xbounds[0], _ybounds[0], _xbounds[1], _ybounds[1], frame_labels.c_str());
            frame->GetXaxis()->CenterTitle(true);
            frame->GetYaxis()->CenterTitle(true);

            if (_custom_z){ graph->SetMinimum(_zbounds[0]); graph->SetMaximum(_zbounds[1]); };
            draw_options += " SAME";
        };

        if (_custom_region)
        {
            // Make an exclusion of the dalitz region
            int Nreg = _region[0].size();
            TCutG *cut = new TCutG("cut", Nreg);
            cut->SetVarX("y");
            cut->SetVarY("x");

            TGraph* outline = new TGraph(Nreg);
            outline->SetLineWidth(2);
            for (int i = 0; i < Nreg; i++) 
            {
                cut->SetPoint    (i, _region[0][i], _region[1][i]);
                outline->SetPoint(i, _region[0][i], _region[1][i]);
            };
            draw_options += " [cut]";

            ex->Draw(); 
            graph->Draw(draw_options.c_str());
            outline->Draw("L SAME");    
            return;
        };

        ex->Draw();
        graph->Draw(draw_options.c_str());
        return;
    };
};