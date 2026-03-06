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
#include "TRandom.h"

using json = nlohmann::json;

namespace iterateKT { namespace COMPASS
{
    // Static identifiers for data_set types
    static const int kReal = 0, kImag = 1, kAbs = 2, kReal1D = 3, kImag1D = 4, kDalitz = 5;
    
    inline std::array<data_set,2>  parse_JSON_1D(std::string input)
    {
        data_set re_out, im_out;
        
        // ---------------------------------------------------------------------------
        // Read in json and organize everything 

        std::string path_to_file = analysis_dir() + "COMPASS_pi1/raw_files/" + input;
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
        double t = -(t_upper + t_lower)/2;
    
        std::string id = "m3π = " + to_string(m3pi,2) + ", t' = " + to_string(-t,2);

        auto bins      = data["bin_centers"];
        auto reF       = data["real(F)"];
        auto imF       = data["imag(F)"];
        int N          = bins.size();

        std::vector<double> s, re, im;
        for (int i = 0; i < N; i++)
        {
            s.push_back(bins[i]);
            re.push_back(reF[i]);
            im.push_back(imF[i]);
        };

        // ---------------------------------------------------------------------------
        //  Organize everything
        re_out._N    = N;                   im_out._N    = N; 
        re_out._id   = id;                  im_out._id   = id;               
        re_out._type = kReal1D;             im_out._type = kImag1D;     
        re_out._extras["Nbins"] = N;        im_out._extras["Nbins"] = N; 
        re_out._extras["m3pi"] = m3pi;      im_out._extras["m3pi"] = m3pi; 
        re_out._extras["t"]    = t;         im_out._extras["t"]    = t;    
        re_out._x = s;                      im_out._x = s;  
        re_out._z = re;                     im_out._z = im;               

        return {re_out, im_out};
    };

    // Parse a JSON file importing everything in a data_set object
    // Columns correspond to: s, t, Abs(M), Err(M)
    inline data_set  parse_JSON(std::string input)
    {
        // Final outputs
        data_set out;

        // ---------------------------------------------------------------------------
        // Read in json and organize everything 

        std::string path_to_file = analysis_dir() + "COMPASS_pi1/raw_files/" + input;
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
        double t = -(t_upper + t_lower)/2;
    
        std::string id = "m3π = " + to_string(m3pi,3) + ", t' = " + to_string(-t,3);

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
                if (are_equal(s1, s2))             continue;
                
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
        out._type = kDalitz;     
        out._extras["Nbins"] = N; 
        out._extras["m3pi"] = m3pi; 
        out._extras["t"]    = t;    
        out._x = sig1;  
        out._y = sig2;             
        out._z = absM; out._dz = errM;               

        return out;
    };

    // Do the above but input bin numbers IDs which are subsequently saved in the data_set
    inline data_set parse_JSON(uint m3pi_bin, uint t_bin)
    {
        std::string sm3pi = to_string(m3pi_bin), st = to_string(t_bin);
        std::string filename = "tBin_"+st+"/dalitz_m3piBin_"+sm3pi+"_tBin_"+st+".json";
        auto out = parse_JSON(filename);
        out._extras["t_bin"]    = t_bin;
        out._extras["m3pi_bin"] = m3pi_bin;
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

        std::string path_to_file = analysis_dir() + "COMPASS_pi1/raw_files/" + input;
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
        double t = -(t_upper + t_lower)/2;
    
        std::string id = "m3π = " + to_string(m3pi,2) + ", -t = " + to_string(-t,2);

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
        out_real._type = kReal;      out_imag._type = kImag;
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

    // ------------------------------------------------------------------------------------------------------------
    // These methods produce pseudodata from the original data file to run a bootstrap

    // Take in a data file from input and resample the data assuming gaussian errors
    inline data_set generate_pseudodata(std::string input, TRandom * rand)
    {
        // Final outputs
        data_set out;

        // ---------------------------------------------------------------------------
        // Read in json and organize everything 

        std::string path_to_file = analysis_dir() + "COMPASS_pi1/raw_files/" + input;
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
        double t = -(t_upper + t_lower)/2;
    
        std::string id = "m3π = " + to_string(m3pi,3) + ", t' = " + to_string(-t,3);

        auto bins      = data["bin_centers"];
        auto abs_M     = data["abs_M"];
        auto std_abs_M = data["std_abs_M"];
        int N          = bins.size();

        // ---------------------------------------------------------------------------
        // Need to filter out any data outside of the physical kinematic region
    
        kinematics kin = new_kinematics(m3pi, M_PION);
        std::vector<double> sig1, sig2, absM, errM;

        // Random number generator with random seed

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                double s1 = bins[i];
                double s2 = bins[j];
                s1 *= s1; s2 *= s2; // mass squared
    
                if (!kin->in_decay_region(s1, s2)) continue;
                if (are_equal(s1, s2))             continue;
                
                sig1.push_back(s1); sig2.push_back(s2);

                // Here instead of saving abs_M, we resample it
                double mean    = abs_M[i][j];
                double std_dev = std_abs_M[i][j];

                double resampled = rand->Gaus(mean, std_dev);
                absM.push_back(resampled); 
                errM.push_back(std_dev);
            };
        };
        int N_actual = sig1.size();

        // ---------------------------------------------------------------------------
        //  Organize everything
        out._N    = N_actual;         
        out._id   = id;               
        out._type = kDalitz;     
        out._extras["Nbins"] = N; 
        out._extras["m3pi"]  = m3pi; 
        out._extras["t"]     = t;    
        out._x = sig1;  
        out._y = sig2;             
        out._z = absM; out._dz = errM;               

        return out; 
    };

    // Do the above but input bin numbers IDs which are subsequently saved in the data_set
    inline data_set generate_pseudodata(uint m3pi_bin, uint t_bin, TRandom * rand)
    {
        std::string sm3pi = to_string(m3pi_bin), st = to_string(t_bin);
        std::string filename = "tBin_"+st+"/dalitz_m3piBin_"+sm3pi+"_tBin_"+st+".json";
        auto out = generate_pseudodata(filename, rand);
        out._extras["t_bin"]    = t_bin;
        out._extras["m3pi_bin"] = m3pi_bin;
        return out;
    };
}; };

#endif 