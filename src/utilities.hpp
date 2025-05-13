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

namespace iterateKT
{
    // ---------------------------------------------------------------------------
    // Angles
    
    inline double degrees(double radians){ return radians*DEG2RAD; };
    inline double radians(double degrees){ return degrees/DEG2RAD; };

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
    inline T kallen(T x, T y, T z)
    {
        return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
    };

    // If any of them are complex, return complex
    inline complex kallen(complex z, double a, double b) { return kallen<complex>(z, complex(1,0)*a, complex(1,0)*b); };
    inline complex kallen(double a, complex z, double b) { return kallen<complex>(complex(1,0)*a, z, complex(1,0)*b); };
    inline complex kallen(double a, double b, complex z) { return kallen<complex>(complex(1,0)*a, complex(1,0)*b, z); };

    inline int sign(double x){ return (x >= 0) ? +1 : -1; }

    // ---------------------------------------------------------------------------
    // Function for easier comparison of doubles using the EPS value defined above
    // be careful when using this in general purposes since its a fixed-tolerance comparision and not always appropriate

    inline bool are_equal(double a, double b)
    {
        return ( abs(a - b) < EPS );
    }

    inline bool are_equal(double a, double b, double tol)
    {
        return ( abs(a - b) < tol );
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
        return (abs(a) < EPS);
    };

    // Aliases for special cases of the above
    template<typename T>
    inline bool is_zero(T a, double tol)
    {
        return (abs(a) < tol);
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

    template <uint NFIRST=PRINT_SPACING, uint NREST=NFIRST, typename First, typename... Rest>
    inline void print(First first, Rest... rest)
    {
        std::cout << std::boolalpha << std::left << std::setprecision(9);  
        std::cout << std::setw(NFIRST) << first;
        print<NREST>(rest...);
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
    inline void print_to_file(std::string outname, std::array<std::vector<double>,N> data, int precision = 12)
    {
        std::ofstream out;
        out.open(outname);
        int spacing = precision + 10;
        for (uint j = 0; j < data[0].size(); j++)
        {
            out << std::left << std::setprecision(precision);
            for (uint i = 0; i < N; i++)
            {
                out << std::setw(spacing) << data[i][j];
            }
            out << std::endl;
        };

        out.close();
        return;
    };

    template<int N>
    inline void print_to_file(std::string outname, std::array<std::string,N> headers, std::array<std::vector<double>,N> data, int precision = 12)
    {
        std::ofstream out;
        out.open(outname);
        int spacing = precision + 10;
        out << std::left << std::setprecision(precision) << std::setw(spacing) << "#" + headers[0];
        for (int i = 1; i < N; i++)
        {
            out << std::setw(spacing) << headers[i]; 
        };
        out << std::endl;

        for (uint j = 0; j < data[0].size(); j++)
        {
            out << std::left;
            for (uint i = 0; i < N; i++)
            {
                out << std::setw(spacing) << data[i][j];
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

    // Another directory that comes up often is the data folder
    inline std::string data_dir(){ return main_dir() + "/data/"; };

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
    template<typename T>
    inline std::vector<T> operator*( std::vector<T> lhs, double c)
    {
        std::vector<T> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i]*c );
        };
        return result;
    };

    template<typename T>
    inline std::vector<T> operator*(double c, std::vector<T> rhs)
    {
        std::vector<T> result;
        for (int i = 0; i < rhs.size(); i++)
        {
            result.push_back( c*rhs[i] );
        };
        return result;
    };

    template<typename T>
    inline std::vector<T> operator/( std::vector<T> lhs, double c)
    {
        std::vector<T> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i]/c );
        };
        return result;
    };

    template<typename T>
    inline std::vector<T> operator-(const std::vector<T> & x)
    {
        return -1. * x;
    };

    template<typename T>
    inline std::vector<T> operator/=(const std::vector<T> & x, double c)
    {
        return x/c;
    };

    template<typename T>
    inline std::vector<T> operator*=(const std::vector<T> & x, double c)
    {
        return x*c;
    };
    
    template<typename T>
    inline std::vector<T> operator+(std::vector<T> lhs, std::vector<T> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<T> result;
        for (uint i = 0; i < lhs.size(); i++) result.push_back( lhs[i] + rhs[i] );
        return result;
    };

    template<typename T>
    inline std::vector<T> operator-(std::vector<T> lhs, std::vector<T> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<double> result;
        for (uint i = 0; i < lhs.size(); i++) result.push_back( lhs[i] - rhs[i] );
        return result;
    };

    // Add a constant to all elements of a vector
    template<typename T>
    inline std::vector<T> operator+(double lhs, std::vector<T> rhs)
    {
        std::vector<T> result;
        for (uint i = 0; i < rhs.size(); i++) result.push_back( lhs + rhs[i] );
        return result;
    };   
    
    template<typename T>
    inline std::vector<T> operator-(double lhs, std::vector<T> rhs) 
    {
        return lhs + (-rhs);
    };

    template<typename T>
    inline std::vector<T> operator+(std::vector<T> lhs, double rhs)
    {
        std::vector<T> result;
        for (uint i = 0; i < lhs.size(); i++) result.push_back( rhs + lhs[i] );
        return result;
    };

    template<typename T>
    inline std::vector<T> operator-(std::vector<T> lhs, double rhs) 
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
    inline std::vector<double> square_elementwise(std::vector<double> in){ return multiply_elementwise(in, in); };

    inline std::vector<double> real(std::vector<complex> vx)
    {
        std::vector<double> out;
        for (auto x : vx) out.push_back(real(x));
        return out;
    };
    inline std::vector<double> imag(std::vector<complex> vx)
    {
        std::vector<double> out;
        for (auto x : vx) out.push_back(real(x));
        return out;
    };

    // ---------------------------------------------------------------------------
    // Simple functions to calculate numerical derivatives up to error ~ h^4
    // Coefficients taken from https://en.wikipedia.org/wiki/Finite_difference_coefficient
    
    // Central difference 
    template<typename T> 
    inline T central_difference_derivative(uint n, std::function<T(double)> F, double x, double h = 1E-3)
    {
        T f3m, f2m, fm, f, fp, f2p, f3p;
        f   = F(x);
        fp  = F(x+h),    fm = F(x-h);
        f2p = F(x+2*h), f2m = F(x-2*h);
        if (n >= 3) { f3p = F(x+3*h); f3m = F(x-3*h); };
        
        T num;
        switch (n)
        {
            case 0  : return F(x);
            case 1  : num =          +f2m/12  -2*fm/3          +2*fp/3   -f2p/12;          break;
            case 2  : num =          -f2m/12  +4*fm/3  -5*f/2  +4*fp/3   -f2p/12;          break;
            case 3  : num = +f3m/8   -f2m    +13*fm/8         -13*fp/8   +f2p     -f3p/8;  break;
            case 4  : num = -f3m/6 +2*f2m    -13*fm/2 +28*f/3 -13*fp/2 +2*f2p     -f3p/6;  break;
            default : 
            {
                warning("central_difference_derivative", "Order n = "+to_string(n)+" derivatives not implemented!"); 
                return NaN<T>();
            };
        };
        return num / pow(h, n);
    };

    // Mixed central derivative of a function of 2 variables d2F(x,y)/dxdy
    template<typename T>
    inline T mixed_partial_derivatives(std::function<T(double,double)> F, std::array<double,2> xs, double e)
    {
        double x = xs[0], y = xs[1];
        T f2p2p = F(x+2*e,y+2*e), f2p2m = F(x+2*e,y-2*e), f2m2p = F(x-2*e,y+2*e), f2m2m = F(x-2*e,y-2*e);
        T f2pp  = F(x+2*e,y+e),   f2pm  = F(x+2*e,y-e),   f2mp  = F(x-2*e,y+e),   f2mm  = F(x-2*e,y-e);
        T fp2p  = F(x+e,y+2*e),   fm2p  = F(x-e,y+2*e),   fp2m  = F(x+e,y-2*e),   fm2m  = F(x-e,y-2*e);
        T fpp   = F(x+e,y+e),     fpm   = F(x+e,y-e),     fmp   = F(x-e,y+e),     fmm   = F(x-e,y-e);
        return (8*(fp2m+f2pm+f2mp+fm2p)     -  8*(fm2m+f2mm+fp2p+f2pp)
                 -(f2p2m+f2m2p-f2m2m-f2p2p) + 64*(fmm+fpp-fpm-fmp))/144/e/e;
    };

    // Forward difference
    template<typename T> 
    inline T forward_difference_derivative(uint n, std::function<T(double)> F, double x, double h = 1E-3)
    {
        if (n == 0) return F(x);

        std::vector<double> c;
        switch (n)
        {
            // O(h^4)
            case 1  : c = {-25./12, 4., -3., 4./3, -1./4}; break;
            case 2  : c = {15./4, -77./6, 107./6, -13., 61./12, -5./6}; break;
            case 3  : c = {-49./8, 29., -461./8, 62., -307./8, 13., -15./8}; break;
            case 4  : c = {28./3, -111./2, 142., -1219/6., 176., -185./2, 82./3 ,-7./2}; break;

            // // O(h^6)
            // case  1 : c = {-49./20, 6., -15./2, 20./3, -15./4, 6./5, -1./6}; break;
            // case  2 : c = {469./90, -223./10, 879./20, -949./18, 41., -201./10, 1019./180, -7./10}; break;
            // case  3 : c = {-801./80, 349./6, -18353./120, 2391./10, -1457./6, 4891./30, -561./8, 527./30, -469./240}; break;
            default : 
            {
                warning("forward_difference_derivative", "Order n = "+to_string(n)+" derivatives not implemented!"); 
                return NaN<T>();
            };
        };
        T num = 0;
        for (int i = 0; i < c.size(); i++) num += c[i] * F(x+i*h);

        return num / pow(h, n);
    };
    
    // and finally backward finite difference
    template<typename T> 
    inline T backward_difference_derivative(uint n, std::function<T(double)> F, double x, double h = 1E-3)
    {
        if (n == 0) return F(x);

        std::vector<double> c;
        switch (n)
        {
            // These are O(h^4)
            case 1  : c = {-25./12, 4., -3., 4./3, -1./4}; break;
            case 2  : c = {15./4, -77./6, 107./6, -13., 61./12, -5./6}; break;
            case 3  : c = {-49./8, 29., -461./8, 62., -307./8, 13., -15./8}; break;
            case 4  : c = {28./3, -111./2, 142., -1219/6., 176., -185./2, 82./3 ,-7./2}; break;

            // // These are O(h^6)
            // case  1 : c = {-49./20, 6., -15./2, 20./3, -15./4, 6./5, -1./6}; break;
            // case  2 : c = {469./90, -223./10, 879./20, -949./18, 41., -201./10, 1019./180, -7./10}; break;
            // case  3 : c = {-801./80, 349./6, -18353./120, 2391./10, -1457./6, 4891./30, -561./8, 527./30, -469./240}; break;
            default : 
            {
                warning("forward_difference_derivative", "Order n = "+to_string(n)+" derivatives not implemented!"); 
                return NaN<T>();
            };
        };
        T num = 0;
        for (int i = 0; i < c.size(); i++) num += c[i] * F(x-i*h);

        return pow(-1, n) * num / pow(h, n);
    };

};
// ---------------------------------------------------------------------------

#endif