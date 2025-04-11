// Methods to interface with COMPASS data sets in JSON format
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef COMPASS_DATA_HPP
#define COMPASS_DATA_HPP

#include "nlohmann/json.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <tuple>
#include "constants.hpp"
#include "kinematics.hpp"
#include "utilities.hpp"
#include "data_set.hpp"

using json = nlohmann::json;

namespace iterateKT { namespace COMPASS
{

    struct fit
    {
        // Static identifiers for data_set types
        static const int kReal = 0, kImag = 1, kAbs = 2;
        static std::string data_type(int i)
        {
            switch (i)
            {
                case kReal: return "Re (M)";
                case kImag: return "Im (M)";
                case kAbs:  return "Abs (M)";
                default: return "ERROR!";
            };
        };

        // Function to minimize
        // Filters whether we're looking at the real or imaginary parts 
        static double fcn(std::vector<data_set> & data_vector, amplitude to_fit)
        {
            double chi2 = 0;
            for (auto data : data_vector)
            {
                for (int i = 0; i < data._N; i++)
                {
                    double from_data  = data._z[i];
                    iterateKT::complex from_model = to_fit->evaluate(data._x[i], data._y[i]);  

                    switch (data._type)
                    {
                        // These two use difference of squares
                        case kReal: chi2 += norm(from_data - real(from_model));        break;
                        case kImag: chi2 += norm(from_data - imag(from_model));        break;
                        // This is a true chi2
                        case kAbs:
                        {
                            if (is_zero(data._dz[i])) continue;
                            chi2  += norm((from_data - abs(from_model)) / data._dz[i]); 
                            break;
                        };
                        default: break;
                    };
                };
            };
            return chi2;
        };
    };

    // Parse a JSON file importing everything in a data_set object
    // Columns correspond to: s, t, Abs(M), Err(M)
    inline data_set  parse_JSON(std::string input)
    {
        // Final outputs
        data_set out;

        // ---------------------------------------------------------------------------
        // Read in json and organize everything 

        std::string path_to_file = data_dir() + "COMPASS/raw_files/" + input;
        std::ifstream raw_file(path_to_file);
        if (!raw_file) fatal("Could not open file: " + path_to_file);
        json data = json::parse(raw_file);
    
        // Calculate central m3pi in bin
        auto m3pi_upper = data["bin_ranges"]["m3pi_upper_limit"].template get<double>();
        auto m3pi_lower = data["bin_ranges"]["m3pi_lower_limit"].template get<double>();
        double m3pi = (m3pi_upper + m3pi_lower)/2;
        
        // Calculate central t in bin
        auto t_upper = data["bin_ranges"]["t_upper_limit"].template get<double>();
        auto t_lower = data["bin_ranges"]["t_lower_limit"].template get<double>();
        double t = (t_upper + t_lower)/2;
    
        std::string id = "m3π = " + to_string(m3pi,2) + ", t' = " + to_string(t,2);

        auto bins      = data["bin_centers"];
        auto abs_M     = data["abs_M"];
        auto std_abs_M = data["std_abs_M"];
        int N          = bins.size();

        // ---------------------------------------------------------------------------
        // Need to filter out any data outside of the physical kinematic region
    
        kinematics kin = new_kinematics(m3pi, M_PION);
        std::vector<double> sig1, sig2, absM, errM;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                double s1 = bins[i];
                double s2 = bins[j];
                s1 *= s1; s2 *= s2; // mass squared
    
                if (!kin->in_decay_region(s1, s2)) continue;
                
                sig1.push_back(s1); sig2.push_back(s2);
                absM.push_back(     abs_M[i][j] ); 
                errM.push_back( std_abs_M[i][j] );
            };
        };
        int N_actual = sig1.size();

        // ---------------------------------------------------------------------------
        //  Organize everything
        out._N    = N_actual;         
        out._id   = id;               
        out._type = fit::kAbs;     
        out._extras["Nbins"] = N; 
        out._extras["m3pi"] = m3pi; 
        out._extras["t"]    = t;    
        out._x = sig1;  
        out._y = sig2;             
        out._z = absM; out._dz = errM;               

        return out;
    };

    // Parse a JSON file importing everything in data_set objects
    // Columns correspond to: s, t, Re(A), Im(A)
    inline std::array<data_set,2>  parse_JSON_ReIm(std::string input)
    {
        // Final outputs
        data_set out_real, out_imag;

        // ---------------------------------------------------------------------------
        // Read in json and organize everything 

        std::string path_to_file = data_dir() + "COMPASS/raw_files/" + input;
        std::ifstream raw_file(path_to_file);
        if (!raw_file) fatal("Could not open file: " + path_to_file);
        json data = json::parse(raw_file);
    
        // Calculate central m3pi in bin
        auto m3pi_upper = data["bin_ranges"]["m3pi_upper_limit"].template get<double>();
        auto m3pi_lower = data["bin_ranges"]["m3pi_lower_limit"].template get<double>();
        double m3pi = (m3pi_upper + m3pi_lower)/2;
        
        // Calculate central t in bin
        auto t_upper = data["bin_ranges"]["t_upper_limit"].template get<double>();
        auto t_lower = data["bin_ranges"]["t_lower_limit"].template get<double>();
        double t = (t_upper + t_lower)/2;
    
        std::string id = "m3π = " + to_string(m3pi,2) + ", -t = " + to_string(t,2);

        auto bin_centers = data["bin_centers"];
        auto real_parts  = data["real(M)"];
        auto imag_parts  = data["imag(M)"];
        int N = bin_centers.size();

        // ---------------------------------------------------------------------------
        // Need to filter out any data outside of the physical kinematic region
    
        kinematics kin = new_kinematics(m3pi, M_PION);
        std::vector<double> sig1, sig2, re, im;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                double s1 = bin_centers[i];
                double s2 = bin_centers[j];
                s1 *= s1; s2 *= s2; // mass squared
    
                if (!kin->in_decay_region(s1, s2)) continue;
                
                double x = real_parts[i][j], y = imag_parts[i][j];
                
                sig1.push_back(s1); sig2.push_back(s2);
                re.push_back(x); im.push_back(y);
            };
        };
        int N_actual = sig1.size();

        // ---------------------------------------------------------------------------
        // Import everything into the data_sets
        out_real._N  = N_actual;          out_imag._N = N_actual;
        out_real._id = id;                out_imag._id = id; 
        out_real._type = fit::kReal;      out_imag._type = fit::kImag;
        out_real._extras["Nbins"] = N;    out_imag._extras["Nbins"] = N; 
        out_real._extras["m3pi"]  = m3pi; out_imag._extras["m3pi"]  = m3pi; 
        out_real._extras["t"]     = t;    out_imag._extras["t"]     = t; 
        out_real._x = sig1;               out_imag._x = sig1;
        out_real._y = sig2;               out_imag._y = sig2;
        out_real._z = re;                 out_imag._z = im;           

        return {out_real, out_imag};
    };

    // Parse a JSON file and output it as a four column ascii file
    // Columns correspond to: s, t, Re(A), Im(A)
    inline void convert_JSON_ReIm_to_ascii(std::string input, std::string output = "")
    {
        auto data = parse_JSON_ReIm(input);
        std::string out = (output == "") ? "dalitz_m3pi_" + to_string(data[0]._extras["m3pi"],2) + "_t_" + to_string(data[0]._extras["t"],2) + ".dat"
                                        : output;
        print_to_file<4>(out, {data[0]._x, data[0]._y, data[0]._z, data[1]._z});
    };

    // Parse a JSON file and output it as a four column ascii file
    // Columns correspond to: s, t, Abs(M), Err(M)
    inline void convert_JSON_to_ascii(std::string input, std::string output = "")
    {
        auto data = parse_JSON(input);
        std::string out = (output == "") ? "dalitz_m3pi_" + to_string(data._extras["m3pi"],4) + "_t_" + to_string(data._extras["t"],2) + ".dat"
                                         : output;
        print_to_file<4>(out, {data._x, data._y, data._z, data._dz});
    };
}; };

#endif 