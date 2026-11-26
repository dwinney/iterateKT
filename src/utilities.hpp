// Useful functions and methods
//
// ---------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cmath>
#include <ios>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <functional>

#include "constants.hpp"
#include "jpacUtils/print.hpp"
#include "jpacUtils/math.hpp"

using namespace jpacUtils;

namespace iterateKT
{
    using jpacUtils::complex;
    using jpacUtils::to_string;

    // ---------------------------------------------------------------------------
    // cin's for operating options

    template<typename T>
    inline void ask_for_value(T & variable, std::string message)
    {
        std::cout << message + " ";
        std::string input = "";
        std::getline( std::cin, input );
        if ( !input.empty() )   
        {
            std::istringstream stream(input);
            stream >> variable;
        };
    };

    template<typename T>
    inline void ask_for_value(T & variable, T default_value, std::string message)
    {
        std::cout << message + " ";
        std::string input = "";
        std::getline( std::cin, input );
        if ( !input.empty() )   
        {
            std::istringstream stream (input);
            stream >> variable;
        }
        else std::cout << default_value << " (default)" << std::endl;
    };

    // ---------------------------------------------------------------------------
    
    // Importing data sets we'll need to be able to find the /data/ directory from the 
    // top level one. Thus we need to be able to access the appropriate environment variable
    inline std::string main_dir()
    {
       // Find the correct data file using the top level repo directory
        std::string top_dir;
        char const * env = std::getenv("ITERATEKT");
        if ( env == NULL || std::string(env) == "" )
        {
            return error("main_dir(): Cannot find environment variable ITERATEKT!", "");
        }
        return std::string(env);  
    };

    // Another directory that comes up often is the data folder
    inline std::string analysis_dir(){ return main_dir() + "/analysis/"; };

    // Same as above but looks for DESKTOP
    inline std::string desktop()
    {
       // Find the correct data file using the top level repo directory
        std::string top_dir;
        char const * env = std::getenv("DESKTOP");
        if ( env == NULL || std::string(env) == "" )
        {
            return error("desktop(): Cannot find environment variable DESKTOP!", "");
        }
        return std::string(env);  
    };
};
// ---------------------------------------------------------------------------

#endif