// Class and methods for handling data sets used for fitting
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef DATA_SET_HPP
#define DATA_SET_HPP

#include "utilities.hpp"

#include <fstream>
#include <sstream>
#include <map>

namespace iterateKT
{
    // Forward declare the option type
    enum class option : unsigned int;

    struct data_set
    {
        // Number of data points
        uint _N = 0;

        // String id 
        std::string _id = "data_set";

        // What is being represented by data
        uint _type;
        
        // vectors to store independent and dependent variables and their errors
        std::vector<double> _x, _y, _z, _dx, _dy, _dz;

        // Store a map in which to store any additional information related to the data_set
        std::map<std::string, double> _extras;

        option _option;

        // If we want a data entry in the legend
        bool _add_to_legend = true;
    };

    // ------------------------------------------------------------------------------
    // Rest are methods to get data_sets more easily

    // Import a set of data with N columns with relative path 
    // and full path main_dir/ + rel_path
    template<int N> 
    inline std::array<std::vector<double>,N> import_data(std::string rel_path, bool quitonfail = false)
    {
        // Check if rel_path starts with a / or not
        // if not we add one
        if (rel_path.front() != '/') rel_path = "/" + rel_path;

        // Add the top level dir path to get full file path
        std::array<std::vector<double>, N> result;
        std::string file_path = main_dir() + rel_path;
        std::ifstream infile(file_path);

        if (!infile.is_open())
        {
            if (quitonfail)
            {
               fatal("import_data", "Cannot open file " + file_path + "!");
            }
            return error("import_data: Cannot open file " + file_path + "!", result);
        };

        // Import data!
        std::string line;
        while (std::getline(infile, line))
        {   
            if (line.empty()) continue;        // skips empty lines
            if (line.front() == '#') continue; // Skip comment lines 
            std::istringstream is(line);   

            for (int i = 0; i < N; i++)
            {
                double x;
                is >> x;
                result[i].push_back(x);
            };
        };
            
        return result;
    };

    // Similar to above except that the data is transposed, i.e. rows are the "categories"
    // and the columns are data points. We specify the number of rows in this case
    template<int N>
    inline std::array<std::vector<double>, N> import_transposed(std::string rel_path)
    {
        // Check if rel_path starts with a / or not
        // if not we add one
        if (rel_path.front() != '/') rel_path = "/" + rel_path;

        // Add the top level dir path to get full file path
        std::array<std::vector<double>, N> result;
        std::string file_path = main_dir() + rel_path;
        std::ifstream infile(file_path);

        if (!infile.is_open())
        {
            return error("import_data: Cannot open file " + file_path + "!", result);
        };

        // Import data!
        for (int i = 0; i < N; i++)
        {   
            std::string line;
            std::getline(infile, line);
            if (line.empty()) continue;        // skips empty lines
            if (line.front() == '#') continue; // Skip comment lines 
            std::istringstream is(line);   

            double x;
            while(is >> x)
            {
                result[i].push_back(x);
            };
        };
            
        return result;
    };

    // If data file has more columns than are actually needed,
    // import with import_data and use this to throw out all but the desired columns
    template<int Nin, int Nout> 
    inline std::array<std::vector<double>,Nout> reshape_data(std::array<std::vector<double>,Nin> data, std::array<int,Nout> to_keep)
    {
        std::array<std::vector<double>, Nout> result;

        for (int i = 0; i < Nout; i++)
        {
            result[i] = data[ to_keep[i] ];
        };
        return result;
    };

    // Make sure all the vectors are the correct size
    template<int S>
    inline int check(std::array<std::vector<double>,S> data, std::string id)
    {
        // Grab the size of the first entry
        int N = data[0].size();
        
        // And compare to the rest
        for (auto column : data)
        {
            if (column.size() != N)
            {
                warning("data_set", "Input vectors of " + id + " have mismatching sizes!");
                return 0;
            };
        };

        return N;
    };
};

#endif