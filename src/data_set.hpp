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

namespace iterateKT
{
    struct data_set
    {
        // Number of data points
        int _N = 0;

        // String id 
        std::string _id = "data_set";

        // What is being represented by data
        int _type;
        
        // vectors to store independent and dependent variables and their errors
        std::vector<double> _x, _y, _z, _dx, _dy, _dz;

        // If we want a data entry in the legend
        bool _add_to_legend = true;
    };

    // ------------------------------------------------------------------------------
    // Rest are methods to get data_sets more easily

    // Import a set of data with N columns with relative path 
    // and full path main_dir/ + rel_path
    template<int N> 
    inline std::array<std::vector<double>,N> import_data(std::string rel_path)
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

    // ---------------------------------------------------------------------------
    // Element-wise operations on data vectors

    // Given two vector<double>s of the same size, calculate the average element wise
    inline std::vector<double> operator*( std::vector<double> lhs, double c)
    {
        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i]*c );
        };
        return result;
    };

    inline std::vector<double> operator*(double c, std::vector<double> rhs)
    {
        std::vector<double> result;
        for (int i = 0; i < rhs.size(); i++)
        {
            result.push_back( c*rhs[i] );
        };
        return result;
    };

    inline std::vector<double> operator/( std::vector<double> lhs, double c)
    {
        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i]/c );
        };
        return result;
    };

    inline std::vector<double> operator-(const std::vector<double> & x)
    {
        return -1 * x;
    };

    
    inline std::vector<double> operator+(std::vector<double> lhs, std::vector<double> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i] + rhs[i] );
        };
        return result;
    };

    inline std::vector<double> operator-(std::vector<double> lhs, std::vector<double> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i] - rhs[i] );
        };
        return result;
    };

    // Add a constant to all elements of a vector
    inline std::vector<double> operator+(double lhs, std::vector<double> rhs)
    {
        std::vector<double> result;
        for (int i = 0; i < rhs.size(); i++)
        {
            result.push_back( lhs + rhs[i] );
        };
        return result;
    };   
    inline std::vector<double> operator-(double lhs, std::vector<double> rhs) 
    {
        return lhs + (-rhs);
    };

    inline std::vector<double> operator+(std::vector<double> lhs, double rhs)
    {
        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( rhs + lhs[i] );
        };
        return result;
    };
    inline std::vector<double> operator-(std::vector<double> lhs, double rhs) 
    {
        return lhs + (-rhs);
    };

    inline std::vector<double> multiply_elementwise(std::vector<double> in1, std::vector<double> in2)
    {
        if (in1.size() != in2.size()) warning("multiply_elementwise()", "Input vectors not the same size!");

        std::vector<double> out;
        for (int i = 0; i < in1.size(); i++)
        {
            out.push_back(in1[i]*in2[i]);
        }
        return out;
    };
    inline std::vector<double> square_elementwise(std::vector<double> in){ return multiply_elementwise(in, in); };
};

#endif