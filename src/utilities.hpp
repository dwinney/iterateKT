// Header file with global phyiscal constants.
// Everything is in GeV unless explicitly stated otherwise.
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
#include <complex>
#include <limits>
#include <ios>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

namespace iteratedOKT
{
    using complex = std::complex<double>;

    // ---------------------------------------------------------------------------
    // Mathematical constants 

    #ifndef PI
        const double PI   = M_PI;
    #endif
    const double DEG2RAD  = (M_PI / 180.);
    const double EPS      = 1.e-7;

    // Unit complex numbers
    const complex XR  (1., 0.);
    const complex I   (0., 1.);
    const complex IEPS(0., EPS);

    // PDG Meson masses in GeV
    const double M_PION      = 0.13957000;
    const double M_KAON      = 0.49367700;
    const double M_ETA       = 0.54753;
    const double M_RHO       = 0.77526;
    const double M_OMEGA     = 0.78265;
    const double M_PHI       = 1.01956;

    // ------------------------------------------------------------------------------
    // // NaN's, 0, and 1 for throwing errors with custom data types

    template<typename T>
    inline T NaN()
    {
        return std::numeric_limits<T>::quiet_NaN();
    }

    template<>
    inline complex NaN()
    {
        return complex(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    }

    template<typename T>
    T zero();

    template<>
    inline complex zero() { return 0; };

    template<typename T> 
    T identity();

    template<>
    inline complex identity() { return 1; };

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
        return c + XR*i;
    };

    inline complex operator+(const int& i, const complex & c)
    {
        return XR*i + c;
    };

    inline complex operator-(const complex&c, const int& i)
    {
        return c - XR*i;
    };

    inline complex operator-(const int& i, const complex & c)
    {
        return XR*i - c;
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
    complex csqrt(T x){ return sqrt(x * XR); };

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
    inline complex Kallen(complex z, double a, double b) { return Kallen<complex>(z, XR*a, XR*b); };
    inline complex Kallen(double a, complex z, double b) { return Kallen<complex>(XR*a, z, XR*b); };
    inline complex Kallen(double a, double b, complex z) { return Kallen<complex>(XR*a, XR*b, z); };

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
    inline bool is_zero(double a)
    {
        return (std::abs(a) < EPS);
    };

    // Aliases for special cases of the above
    inline bool is_zero(double a, double tol)
    {
        return (std::abs(a) < tol);
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
    const std::string UNIT_DIV = std::string(PRINT_SPACING, '-');

    // ---------------------------------------------------------------------------   
    // Output an empty line to the terminal
    inline void line()
    {
        std::cout << std::endl;
    };

    // Print out a horizontal line
    inline void divider()
    {
        std::cout << std::string(TEXT_WIDTH, '-') << std::endl;
    };

    inline void divider(int n)
    {
        std::string div;
        for (int i = 0; i < n; i++)
        {
            div = div + UNIT_DIV;
        }
        std::cout << div << std::endl;
    };
    
    inline void dashed_divider()
    {
        std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
    };

    template<typename T>
    inline void print(T x)
    {
        std::cout << std::boolalpha << std::left << std::setprecision(9);  
        std::cout << std::setw(PRINT_SPACING) << x << std::endl;
    };

    template <typename First, typename... Rest>
    inline void print(First first, Rest... rest)
    {
        std::cout << std::boolalpha << std::left << std::setprecision(9);  
        std::cout << std::setw(PRINT_SPACING) << first;
        print(rest...);
    } 

    template<typename T>
    inline void print(std::vector<T> v)
    {
        std::cout << std::boolalpha << std::setprecision(9);  
        for (auto vi : v)
        {
            std::cout << std::left << std::setw(PRINT_SPACING) << vi << std::endl;
        };
    };

    // ---------------------------------------------------------------------------
    // String operations

    // Produce a string with the format "name = value units"

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

        for (int j = 0; j < data[0].size(); j++)
        {
            out << std::left;
            for (int i = 0; i < N; i++)
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

        for (int j = 0; j < data[0].size(); j++)
        {
            out << std::left;
            for (int i = 0; i < N; i++)
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
        char const * env = std::getenv("iteratedOKT");
        if ( env == NULL || std::string(env) == "" )
        {
            return error("main_dir(): Cannot find environment variable iteratedOKT!", "");
        }
        return std::string(env);  
    };

};
// ---------------------------------------------------------------------------

#endif