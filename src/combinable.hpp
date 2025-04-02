// Abstract class to interface with the plotter::combine method
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef COMBINABLE_HPP
#define COMBINABLE_HPP

namespace iterateKT
{
    class combinable
    {
        public:

        combinable(){};

        virtual void combine_draw(double scale) = 0;
    };
};

#endif