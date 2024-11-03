// Simple class to time processes
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include "utilities.hpp"

namespace iterateKT
{
    class timer
    {
        public: 

        timer(bool messages = true) : _messages(messages)
        {};

        ~timer()
        {
            if (_started) warning("timer", "Started timer going out of scope!");
        };

        inline void start()  
        { 
            if (_messages) print("Timer started!");
            _start = std::chrono::high_resolution_clock::now(); 
            _lap   = _start;
            _started = true;
        };
        inline void stop()   
        { 
            if (_messages) print("Timer stopped!");
            _end   = std::chrono::high_resolution_clock::now(); 
            _started = false;
        };
        inline auto lap(std::string message = "")
        {
            auto now = std::chrono::high_resolution_clock::now();
            auto x = std::chrono::duration_cast< std::chrono::seconds>(now - _lap).count();
            std::string extra = (message == "") ? "" : " (" + message + ")";
            if (_messages) print("Lap time" + extra + ": " + std::to_string(x) + "s!");
            _lap = now;
            return x;
        };
        inline auto elapsed()
        {
            auto x = std::chrono::duration_cast< std::chrono::seconds>(_end - _start).count(); 
            return x;
        };
            
        inline void print_elapsed(){ print("Time elapsed: " + std::to_string(elapsed()) + "s!"); };

        private: 
        
        bool _messages = true;
        bool _started  = false;
        std::chrono::time_point<std::chrono::high_resolution_clock> _start, _end, _lap;



    }; //timer
}; //namespace

#endif // TIMER_HPP