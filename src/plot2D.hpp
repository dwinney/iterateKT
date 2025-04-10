// Reworking of the jpacStyle as internal library special for jpacPhoto
// This defines a single plot object, which is what actually prints out the files
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef PLOT2D_HPP
#define PLOT2D_HPP

#include <array>
#include <vector>
#include <deque>
#include <iostream>   
#include <sstream> 
#include <functional>

#include <TROOT.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraph2D.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TExec.h>

#include "data_set.hpp"
#include "colors.hpp"

namespace iterateKT
{
    class plotter;

    // This class contains the entries, data, and options of producing a single plot/file
    // These can be generated from the plotter->make_plot() method which applies
    // all global settings
    class plot2D
    {
        public:

        // Outputting function, generates plot and saves it to file
        void save(); 
        inline void save(std::string filename){ _filename = filename; save(); }; 

        // make sure _data is emptied out
        inline void clear_data(){ for (auto dim : _data) dim.clear(); };

        // Set and clear data
        inline void set_data(std::array<std::vector<double>,3> data){ clear_data(); _data = data; };
        inline void set_data(data_set data){ clear_data(); set_data({data._x, data._y, data._z}); };

        
        // Custom plotting region
        inline void set_region(std::array<std::vector<double>,2> region){ _custom_region = true; _region = region; };
        
        // Set title and axis labels
        inline void set_title(std::string title){ _title = title; };
        inline void set_labels(std::string xl, std::string yl){ _xlabel = xl; _ylabel = yl; };
        
        // Set a custom plotting range
        inline void set_ranges(std::array<double,2> xr, std::array<double,2> yr){ _custom_ranges = true; _xbounds = xr; _ybounds = yr; };
        inline void set_ranges(std::array<double,2> xr, std::array<double,2> yr, std::array<double,2> zr)
        {
            _custom_ranges = true; _xbounds = xr; _ybounds = yr; 
            _custom_z = true;      _zbounds = zr; 
        };
        
        inline void set_Nbins(int x, int y = 100){ _nbins = x; _ncontours = y; };
        inline void set_palette(int x, bool invert = false){ _palette = x; _inverted = invert; };

        private: 

        // Constructor is private, only creatable through plotter
        plot2D(TCanvas* canvas, std::string filename)
        : _canvas(canvas), _filename(filename)
        {};
        friend class plotter;

        // Default filename
        std::string _filename = "";

        // Canvas that the plot actually gets drawn on
        TCanvas* _canvas;
        
        // Save the data to plot
        std::array<std::vector<double>,3> _data;

        // Whether to trim the square plotting region
        bool _custom_region = false;
        std::array<std::vector<double>,2> _region;
        
        // Plot title and axis labels
        std::string _title = "";
        std::string _xlabel = "", _ylabel = "";

        // Custom ranges for both x and y axis
        bool _custom_ranges = false, _custom_z = false;
        std::array<double,2> _xbounds, _ybounds, _zbounds;

        // Related to color options
        int _palette   = kBird; 
        bool _inverted = false;

        int _nbins = 300, _ncontours = 256;
        
        // Apply all settings and draw onto the saved _canvas
        void draw();
    };
};

#endif