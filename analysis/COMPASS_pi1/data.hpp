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
    
    // Parse a JSON file importing everything in a data_set object
    // Columns correspond to: s, t, Abs(M), Err(M)
    inline data_set  parse_JSON(uint m3pibin, std::string input)
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
        std::string bin = "m3pi_bin_number_" + to_string(m3pibin);
        auto m3pi_upper = data["bins"][bin]["bin_ranges"]["m3pi_upper_limit"].template get<double>();
        auto m3pi_lower = data["bins"][bin]["bin_ranges"]["m3pi_lower_limit"].template get<double>();
        double m3pi = (m3pi_upper + m3pi_lower)/2;
        
        // Calculate central t in bin
        auto t_upper = data["bins"][bin]["bin_ranges"]["t_upper_limit"].template get<double>();
        auto t_lower = data["bins"][bin]["bin_ranges"]["t_lower_limit"].template get<double>();
        double t = -(t_upper + t_lower)/2;
    
        std::string id = "m3π = " + to_string(m3pi,3) + ", t' = " + to_string(-t,3);

        auto bins      = data["bins"][bin]["bin_centers"];
        auto abs_M     = data["bins"][bin]["abs_M"];
        auto std_abs_M = data["bins"][bin]["std_abs_M"];
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

                double z = abs_M[i][j];
                if (is_zero(z))                    continue;
                
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
        out._extras["tbin_width"]    = (t_upper - t_lower)/2;
        out._extras["m3pibin_width"] = (m3pi_upper - m3pi_lower)/2;
        out._x = sig1;  
        out._y = sig2;             
        out._z = absM; out._dz = errM;               

        return out;
    };

    // Do the above but input bin numbers IDs which are subsequently saved in the data_set
    inline data_set parse_JSON(uint m3pi_bin, uint t_bin)
    {
        std::string st = to_string(t_bin);
        std::string filename = "tBin_"+st+".json";
        auto out = parse_JSON(m3pi_bin, filename);
        out._extras["t_bin"]    = t_bin;
        out._extras["m3pi_bin"] = m3pi_bin;
        return out;
    };

    // Parse a JSON file importing everything in data_set objects
    // Columns correspond to: s, t, Re(A), Im(A)
    inline std::array<data_set,2>  parse_JSON_ReIm(uint m3pibin, std::string input)
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
        std::string bin = "m3pi_bin_number_" + to_string(m3pibin);
        auto m3pi_upper = data["bins"][bin]["bin_ranges"]["m3pi_upper_limit"].template get<double>();
        auto m3pi_lower = data["bins"][bin]["bin_ranges"]["m3pi_lower_limit"].template get<double>();
        double m3pi = (m3pi_upper + m3pi_lower)/2;
        
        // Calculate central t in bin
        auto t_upper = data["bins"][bin]["bin_ranges"]["t_upper_limit"].template get<double>();
        auto t_lower = data["bins"][bin]["bin_ranges"]["t_lower_limit"].template get<double>();
        double t = -(t_upper + t_lower)/2;
    
        std::string id = "m3π = " + to_string(m3pi,2) + ", -t = " + to_string(-t,2);

        auto bin_centers = data["bins"][bin]["bin_centers"];
        auto real_parts  = data["bins"][bin]["real(M)"];
        auto imag_parts  = data["bins"][bin]["imag(M)"];
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
}; };

#endif 