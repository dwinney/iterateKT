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

#include "constants.hpp"

namespace iterateKT
{
    // ---------------------------------------------------------------------------
    // Out of the box std::complex<double> doesnt play well with ints. Here we explicitly
    // give it the functionality we want

    // Additionally the complex type is a liitle dim in c++ so we need to define int & bool multiplication
    inline complex operator*(const int& c, const complex& rhs)
    {
        return complex(c*rhs.real(), c*rhs.imag());
    };

    inline complex operator*(const complex& lhs, const int& c)
    {
        return complex(c*lhs.real(), c*lhs.imag());
    };

    inline complex operator*(const bool& c, const complex& rhs)
    {
        return (c) ? rhs : 0.;
    };

    inline complex operator*(const complex& lhs, const bool& c)
    {
        return (c) ? lhs : 0.;
    };

    inline complex operator/(const complex&c, const int& i)
    {
        return (1./i)*c;
    };

    inline complex operator/(const int& i, const complex&c)
    {
        return (1./c)*i;
    };

    inline complex operator+(const complex&c, const int& i)
    {
        return c + complex(1,0)*i;
    };

    inline complex operator+(const int& i, const complex & c)
    {
        return complex(1,0)*i + c;
    };

    inline complex operator-(const complex&c, const int& i)
    {
        return c - complex(1,0)*i;
    };

    inline complex operator-(const int& i, const complex & c)
    {
        return complex(1,0)*i - c;
    };

    inline bool operator == (const complex &z,const int n)
    {
        return (z == static_cast<double> (n));
    }
    inline bool operator == (const int n, const complex &z)
    {
        return (z == static_cast<double> (n));
    }
    inline bool operator != (const complex &z,const int n)
    {
        return (z != static_cast<double> (n));
    }
    inline bool operator != (const int n, const complex &z)
    {
        return (z != static_cast<double> (n));
    }

    // This makes it so we always default to complex regardless of whether the input is an int or double
    template<typename T>
    complex csqrt(T x){ return sqrt(x * complex(1,0)); };

    inline unsigned int factorial(unsigned int n) 
    {
        if (n == 0)
        return 1;
        return n * factorial(n - 1);
    };

    // ---------------------------------------------------------------------------
    // Kallen Triangle function

    // Only way to get a double or int Kallen is if all inputs are double/int
    template<typename T>
    inline T Kallen(T x, T y, T z)
    {
        return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
    };

    // If any of them are complex, return complex
    inline complex Kallen(complex z, double a, double b) { return Kallen<complex>(z, complex(1,0)*a, complex(1,0)*b); };
    inline complex Kallen(double a, complex z, double b) { return Kallen<complex>(complex(1,0)*a, z, complex(1,0)*b); };
    inline complex Kallen(double a, double b, complex z) { return Kallen<complex>(complex(1,0)*a, complex(1,0)*b, z); };

    inline int sign(double x){ return (x >= 0) ? +1 : -1; }

    // ---------------------------------------------------------------------------
    // Function for easier comparison of doubles using the EPS value defined above
    // be careful when using this in general purposes since its a fixed-tolerance comparision and not always appropriate

    inline bool are_equal(double a, double b)
    {
        return ( std::abs(a - b) < EPS );
    }

    inline bool are_equal(double a, double b, double tol)
    {
        return ( std::abs(a - b) < tol );
    }

    // Same thing for comparing complex doubles
    inline bool are_equal(complex a, complex b)
    {
        return (are_equal(real(a), real(b)) && are_equal(imag(a), imag(b)));
    };

    // Same thing for comparing complex doubles
    inline bool are_equal(complex a, complex b, double tol)
    {
        return (are_equal(real(a), real(b), tol) && are_equal(imag(a), imag(b), tol));
    };

    // Aliases for special cases of the above
    template<typename T>
    inline bool is_zero(T a)
    {
        return (std::abs(a) < EPS);
    };

    // Aliases for special cases of the above
    template<typename T>
    inline bool is_zero(T a, double tol)
    {
        return (std::abs(a) < tol);
    };

    inline bool is_odd(int a)
    {
        return (a%2 == 1);
    };

    // ---------------------------------------------------------------------------
    // ERROR Messages
    
    // Throw an error message then quits code 
    inline void fatal()
    {
        std::cout << std::left << "FATAL ERROR! Quitting..." << std::endl;
        exit( EXIT_FAILURE );
    };

    // Error message with location and reason messages too
    inline void fatal(std::string location, std::string reason = "")
    {
        std::cout << std::left << "FATAL ERROR! " + location + ": " + reason << std::endl;
        std::cout << std::left << "Quitting..." << std::endl;

        exit( EXIT_FAILURE );
    };

    // Warning message does not exit code or returns simply throws a message up
    inline void warning(std::string message)
    {
        std::cout << std::left << "WARNING! " + message << std::endl;
    };

    // Warning message with additional location
    inline void warning(std::string location, std::string message)
    {
        std::cout << std::left << "WARNING! " + location + ": " + message << std::endl;
    };

    // Throw an error message without location and return a value
    template<typename T> 
    inline T error(std::string location, std::string message, T return_value )
    {
        warning(location, message);
        return return_value;
    };
    
    template<typename T> 
    inline T error(std::string message, T return_value )
    {
        warning(message);
        return return_value;
    };

    // Alternatively without a return value, this simply returns void type
    inline void error(std::string message)
    {
        warning(message);
        return;
    };
   
    // ---------------------------------------------------------------------------   
    // Default values
    const int TEXT_WIDTH       = 62;
    const int PRINT_SPACING    = 15;
    const int PRINT_PRECISION  = 9;    
    const int STRING_PRECISION = 3;

    // ---------------------------------------------------------------------------   
    // Output an empty line to the terminal
    inline void line()
    {
        std::cout << std::endl;
    };

    // Print out a horizontal line
    template<uint N=TEXT_WIDTH>
    inline void divider()
    {
        std::cout << std::string(N, '-') << std::endl;
    };

    template<uint N=PRINT_SPACING>
    inline void divider(int n)
    {
        std::string unit_div = std::string(N, '-');
        std::string div;
        for (int i = 0; i < n; i++)
        {
            div = div + unit_div;
        }
        std::cout << div << std::endl;
    };
    
    inline void dashed_divider()
    {
        std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
    };

    template<uint N=PRINT_SPACING, typename T>
    inline void print(T x)
    {
        std::cout << std::boolalpha << std::left << std::setprecision(9);  
        std::cout << std::setw(N) << x << std::endl;
    };

    template <uint N=PRINT_SPACING, typename First, typename... Rest>
    inline void print(First first, Rest... rest)
    {
        std::cout << std::boolalpha << std::left << std::setprecision(9);  
        std::cout << std::setw(N) << first;
        print(rest...);
    } 

    template<uint N=PRINT_SPACING, typename T>
    inline void print(std::vector<T> v)
    {
        std::cout << std::boolalpha << std::setprecision(9);  
        for (auto vi : v)
        {
            std::cout << std::left << std::setw(N) << vi << std::endl;
        };
    };

    // ---------------------------------------------------------------------------
    // String operations

    // Produce a string with the format "name = value units"

    inline std::string to_string(double d, uint precision = 8)
    {
        std::stringstream ss;
        ss << std::setprecision(precision) << d;
        return ss.str();
    };

    template <typename T>
    inline std::string var_def(std::string name, T value, std::string units = "")
    {
        std::stringstream ss;
        ss << std::setprecision(STRING_PRECISION) << name + " = " << value << " " + units;
        return ss.str();
    };

    // Print a string centered on the terminal 
    inline void centered(int n, std::string words)
    {
        int x = words.length();
        int gap_width = (n * PRINT_SPACING - x)/2;
        std::cout << std::left << std::setw(gap_width) << "" << std::setw(x) << words << std::setw(gap_width) << "" << std::endl;
    };

    // ---------------------------------------------------------------------------
    // Print N columns of data to file
    // We assume theyre all the same size, if not this breaks and its not our fault

    template<int N>
    inline void print_to_file(std::string outname, std::array<std::vector<double>,N> data)
    {
        std::ofstream out;
        out.open(outname);

        for (uint j = 0; j < data[0].size(); j++)
        {
            out << std::left;
            for (uint i = 0; i < N; i++)
            {
                out << std::setw(PRINT_SPACING) << data[i][j];
            }
            out << std::endl;
        };

        out.close();
        return;
    };

    template<int N>
    inline void print_to_file(std::string outname, std::array<std::string,N> headers, std::array<std::vector<double>,N> data)
    {
        std::ofstream out;
        out.open(outname);

        out << std::left << std::setw(PRINT_SPACING) << "#" + headers[0];
        for (int i = 1; i < N; i++)
        {
            out << std::setw(PRINT_SPACING) << headers[i]; 
        };
        out << std::endl;

        for (uint j = 0; j < data[0].size(); j++)
        {
            out << std::left;
            for (uint i = 0; i < N; i++)
            {
                out << std::setw(PRINT_SPACING) << data[i][j];
            }
            out << std::endl;
        };

        out.close();
        return;
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

    inline std::vector<double> operator/=(const std::vector<double> & x, double c)
    {
        return x/c;
    };

    inline std::vector<double> operator*=(const std::vector<double> & x, double c)
    {
        return x*c;
    };
    
    inline std::vector<double> operator+(std::vector<double> lhs, std::vector<double> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<double> result;
        for (uint i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i] + rhs[i] );
        };
        return result;
    };

    inline std::vector<double> operator-(std::vector<double> lhs, std::vector<double> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<double> result;
        for (uint i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i] - rhs[i] );
        };
        return result;
    };

    // Add a constant to all elements of a vector
    inline std::vector<double> operator+(double lhs, std::vector<double> rhs)
    {
        std::vector<double> result;
        for (uint i = 0; i < rhs.size(); i++)
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
        for (uint i = 0; i < lhs.size(); i++)
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
        for (uint i = 0; i < in1.size(); i++)
        {
            out.push_back(in1[i]*in2[i]);
        }
        return out;
    };
    inline std::vector<double> square(std::vector<double> in){ return multiply_elementwise(in, in); };

};
// ---------------------------------------------------------------------------

#endif