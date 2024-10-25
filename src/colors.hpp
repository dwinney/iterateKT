// TColor inex definitions for the JPAC color palatte used for plotting
// These are not actually defined here instead 
// they are initialized by the plotter object in plotter.hpp
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef COLORS_HPP
#define COLORS_HPP

#include <TColor.h>

namespace iterateKT
{
    // Give each jpacColor a ROOT TColor index name so it may be called globally after its initialized
    enum class jpacColor: Int_t { Blue   = 2001, Red      = 2002, Green = 2003,
                                  Orange = 2004, Purple   = 2005, Brown = 2006, 
                                  Pink   = 2007, Gold     = 2008, Aqua  = 2009, 
                                  Grey   = 2010, DarkGrey = 2011                };

    // Convert from jpacColor to its underlying int
    inline constexpr Int_t operator+(jpacColor x)
    {
        return static_cast<Int_t>(x);
    };

    constexpr std::array<jpacColor,11> JPACCOLORS = {jpacColor::Blue,   jpacColor::Red,    jpacColor::Green, 
                                                     jpacColor::Orange, jpacColor::Purple, jpacColor::Brown, 
                                                     jpacColor::Pink,   jpacColor::Gold,   jpacColor::Aqua, 
                                                     jpacColor::Grey,   jpacColor::DarkGrey };

                                                     
    struct entry_style
    {
        jpacColor _color      = jpacColor::DarkGrey;  // Color code 
        int  _style           = 0;                    // Either linestyle or markerstyle code
        bool _add_to_legend   = true;                 // Whether to add this curve to the legend
        std::string _label    = "";                   // Label to add to Legend
        std::string _draw_opt = "L";                // string which enters ROOT::Draw() 
    };

    inline entry_style dashed(jpacColor color, std::string id = "")
    {
        entry_style dashed;
        dashed._color = color;
        dashed._style = kDashed;
        dashed._label = id;
        dashed._add_to_legend = (dashed._label != "");
        return dashed;
    };

    inline entry_style solid(jpacColor color, std::string id = "")
    {
        entry_style dashed;
        dashed._color = color;
        dashed._style = kSolid;
        dashed._label = id;
        dashed._add_to_legend = (dashed._label != "");
        return dashed;
    };

    inline entry_style dotted(jpacColor color, std::string id = "")
    {
        entry_style dotted;
        dotted._color = color;
        dotted._style = kDotted;
        dotted._label = id;
        dotted._add_to_legend = (dotted._label != "");
        return dotted;
    };

};


#endif 