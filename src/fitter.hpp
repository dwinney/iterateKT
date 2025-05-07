// Class to fit amplitudes to data
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef FITTER_HPP
#define FITTER_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom.h"

namespace iterateKT
{
    // Each free parameter of a model has associated with is a bunch of data
    class parameter
    {
        public:

        // Default initialized
        // This initializes a fixed parameter for normalizations
        parameter()
        : _mod_fixed(true), _arg_fixed(true), 
          _mod(1.), _arg(0.),
          _label(default_label(0)), _i(-1)
        {};
        
        // Constructor for amplitude variables
        parameter(int i) : _i(i), _label(default_label(i)){};

        uint         _i;
        std::string _label;
        std::string _mod_message = "[ ≥ 0 ]", _arg_message = "[-π, π]";
        bool   _mod_fixed         = false;
        bool   _arg_fixed         = false;
        bool   _real              = false;
        std::string mod_label(){ return "|" + _label + "|";    };
        std::string arg_label(){ return "arg(" + _label + ")"; };

        // Each parameter is complex so store both mod and phase seperately
        double _mod           = 0;
        double _arg           = 0;
        complex value(){ return _mod*exp(I*_arg); };

        // Step-size for both
        double _step          = 0.1;

        static inline std::string default_label(int i){ return "par[" + std::to_string(i) + "]"; };
    };

    // ---------------------------------------------------------------------------
    // Actual fitter object
    // This is templated because it requires an implementation of the chi2 function to be fit
    // The template F should contain the details of the fit and the following static functions
    // fcn(amplitude, std::vector<data_set>&) [function to be minimized, e.g. chi2]
    // 
    template<class F>
    class fitter
    {
        public: 

        // Basic constructor, only requires amplitude to be fit 
        // uses default settings for minuit
        fitter(amplitude amp_to_fit)
        : _amplitude(amp_to_fit),
          _minuit(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined"))
        { reset_parameters(); };

        // Parameterized constructor 
        // with explicit choice of minimization strategy and tolerance of minuit routines
        fitter(amplitude amp_to_fit, std::string strategy, double tolerance = 1.E-6)
        : _amplitude(amp_to_fit), _tolerance(tolerance),
          _minuit(ROOT::Math::Factory::CreateMinimizer("Minuit2", strategy))
        { reset_parameters(); }; 

        // -----------------------------------------------------------------------
        // Methods to add/remove data to be fit against

        inline void clear_data(){ _data.clear(); _N = 0; };
        inline void add_data(data_set data){ _N += data._N; _data.push_back(data); };
        inline void add_data(std::vector<data_set> data){ for (auto datum : data) add_data(datum); };
        template <int N>
        inline void add_data(std::array<data_set,N> data){ for (auto datum : data) add_data(datum); };

        // -----------------------------------------------------------------------
        // Set limits, labels, and fix parameters

        // Reset labels, limits and options on all parameters
        inline void reset_parameters()
        {
            std::vector<std::string> labels; // Default label names
            _pars.clear(); _Nfree = 2*_amplitude->N_pars();
            for (int i = 0; i < _amplitude->N_pars(); i++)
            {
                _pars.push_back(i);
                labels.push_back("par[" + std::to_string(i) + "]");
            };
            set_parameter_labels(labels);
        };

        // Give each parameter a label beyond their default par[i] name
        inline void set_parameter_labels(std::vector<std::string> labels)
        {
            if (labels.size() != _pars.size())
            {
                warning("fitter::set_parameter_labels", "Labels vector does not match number of parameters!");
                return;
            }
            for (int i = 0; i < _pars.size(); i++) _pars[i]._label = labels[i];
        };

        inline void fix_argument(parameter & par, double val)
        {
            if (!par._arg_fixed) _Nfree--;
            par._arg_fixed   = true;
            par._arg         = val;
            par._arg_message = "[FIXED]";
        };

        inline void fix_argument(std::string label, double val)
        {
            int index = find_parameter(label);
            if (index < 0) return;
            return fix_argument(_pars[index], val);
        };

        inline void fix_modulus(parameter & par, double val)
        {
            if (!par._mod_fixed) _Nfree--;
            par._mod_fixed   = true;
            par._mod         = val;
            par._mod_message = "[FIXED]";
        };

        inline void fix_modulus(std::string label, double val)
        {
            int index = find_parameter(label);
            if (index < 0) return;
            return fix_modulus(_pars[index], val);
        };

        inline void fix_parameter(parameter & par, complex val)
        {
            fix_modulus(par, abs(val)); fix_argument(par, arg(val));
        };

        inline void fix_parameter(std::string label, complex val)
        {
            int index = find_parameter(label); 
            if (index < 0) return;
            fix_modulus(_pars[index], abs(val)); fix_argument(_pars[index], arg(val));
        };

        inline void make_real(parameter & par)
        {
            par._real = true; 
            par._mod_message = "[REAL]";
            fix_argument(par, 0.);
        };
        inline void make_real(std::string label)
        {
            int index = find_parameter(label); 
            if (index < 0) return;
            make_real(_pars[index]);
        };

        inline void free_parameter(parameter& par)
        {
            par._real = false;
            // if not fixed, this does nothing
            if (par._mod_fixed)
            {
                par._mod_fixed   = false;
                par._mod_message = "[ ≥ 0 ]";
                _Nfree++;
            };
            if (par._arg_fixed)
            {
                par._arg_fixed   = false;
                par._arg_message = "[-π, π]";
                _Nfree++;
            };
        };

        inline void free_parameter(std::string label)
        {
            int index = find_parameter(label);
            if (index < 0) return;
            free_parameter(_pars[index]);
        };

        // Actually do the fit given a vector of size amp->N_pars() as starting values
        // Prints results to command line but also returns the best-fit chi2 value
        inline void do_fit(std::vector<complex> starting_guess, bool show_data = true)
        {
            int expected_size = 0;
            for (auto par : _pars) if (!par._mod_fixed || !par._arg_fixed) expected_size++;
            if (starting_guess.size() !=  expected_size) 
            {
                warning("fitter::do_fit", "Starting guess not the correct size! Expected " + std::to_string(expected_size) + " parameters!");
                return;
            };

            set_up(starting_guess);

            if (show_data) { line(); data_info(); };
            parameter_info();

            auto start = std::chrono::high_resolution_clock::now();
            std::cout << "Beginning fit..." << std::flush; 

            if (_print_level != 0) line();   
            _minuit->Minimize();
            if (_print_level != 0) line();   

            std::cout << "Done! \n";

            // Timing info
            auto stop     = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast< std::chrono::seconds>(stop - start);
            std::cout << std::left << "Finished in " << duration.count() << " s" << std::endl;

            line();
            print_results();
        };

        // -----------------------------------------------------------------------
        // Methods related to fit options

        // Set the maximum number of calls minuit will do
        inline void set_max_calls(int n){ _max_calls = n; };
        
        // Message level for minuit (0-4)
        inline void set_print_level(int n){ _print_level = n; };

        // Change tolerance
        inline void set_tolerance(double tol){ _tolerance = tol; };

        // Number of degrees of freedom
        inline int dof(){ return _N - _minuit->NFree(); };

        private:

        // This ptr should point to the amplitude to be fit
        amplitude _amplitude = nullptr;

        // -----------------------------------------------------------------------
        // Data handling

        int _N = 0;  // Total number of data points
        std::vector<data_set> _data;  // Contain all data

        // -----------------------------------------------------------------------
        // MINUIT handling 

        int _print_level   = 0;     // Error code for MINUIT
        int _max_calls     = 1E6;   // Max calls allowed for minimization fcn
        double _tolerance  = 1.E-6; // Minimization tolerance

        ROOT::Math::Minimizer * _minuit;
        ROOT::Math::Functor     _wfcn;

        // Random number generator for creating initial guesses;
        TRandom *_guesser = new TRandom(0);

        // Initialize minuit with all our parameter options etc
        inline void set_up(std::vector<complex> starting_guess)
        {
            _minuit->Clear();
            _minuit->SetTolerance(_tolerance);
            _minuit->SetPrintLevel(_print_level);
            _minuit->SetMaxFunctionCalls(_max_calls);

            // Iterate over each _par but also keep track of the index in starting_guess 
            // because parameters might be fixed, these indexes dont necessarily line up
            int i = 0, j = 0;
            for (auto &par : _pars)
            {   
                bool move_up = false;
                if (par._real && !par._mod_fixed)
                {
                    // Set up the starting guess
                    par._mod = real(starting_guess[i]);
                    _minuit->SetVariable(j,   par._label,    real(starting_guess[i]), par._step);
                    j++; move_up = true;
                };
                if (!par._real && !par._mod_fixed)
                {
                    // Set up the starting guess
                    par._mod = abs(starting_guess[i]);
                    _minuit->SetVariable(j,   "|" + par._label + "|",    abs(starting_guess[i]), par._step);
                    _minuit->SetVariableLowerLimit(j, 0.);    // mod is positive definite
                    j++; move_up = true;
                }
                if (!par._arg_fixed)
                {   
                    par._arg = std::arg(starting_guess[i]);
                    _minuit->SetVariable(j, "arg(" + par._label + ")", std::arg(starting_guess[i]), par._step);

                    // We DONT limit the angle becasue the fitter does not do well with the 2pi jumps. 
                    // Instead we fit with continuous angle (within reason), and at the end map it back to [-pi,pi]
                    _minuit->SetVariableLimits(j, -3*PI, 3*PI);
                    j++; move_up = true;
                }
                if (move_up) i++;
            };
        
            _wfcn = ROOT::Math::Functor(this, &fitter::fit_fcn, _Nfree);
            _minuit->SetFunction(_wfcn);
        };

                
        // -----------------------------------------------------------------------
        // Calcualtions of chi-squared      
      
        // This is the actual function that gets called by minuit
        inline double fit_fcn(const double *cpars)
        { 
            // First convert the C string to a C++ vector
            std::vector<complex> pars = complex_convert(cpars);
        
            // Sometimes we want the amplitude to do something to the fitter output
            // Before we actually save them, run through amplitudes processing function
            // By default this does nothing
            std::vector<complex> processed = _amplitude->process_fitter_parameters(pars);
    
            // Pass parameters to the amplitude
            _amplitude->set_parameters(processed);

            // Pass both this and data to fit function
            return F::fcn(_data, _amplitude); 
        };
        
        // -----------------------------------------------------------------------
        // Parameter handling

        // Store of parameter info
        std::vector<parameter> _pars;

        // Number of free parameters
        int _Nfree = 0;

        // Default guess_range to initalize parameters
        std::array<double,2> _guess_range = {0, 5};

        // Given a parameter label, find the corresponding index
        inline int find_parameter(std::string label)
        {
            for (auto par : _pars) if (par._label == label) return par._i;
            return error("fitter::find_parameter", "Cannot find parameter labeled " + label + "!", -1);
        };

        // Given a C-style array of size 2*_Nfree
        // Convert to a C++ style std::vector<double> of size 2*_N_free
        inline std::vector<double> convert(const double * cpars)
        {
            std::vector<double> results;
            for (int i = 0; i < _minuit->NFree(); i++) results.push_back(cpars[i]);
            return results;
        };
        
        // Given a C-style array of size 2*_Nfree
        // Convert to a C++ style std::vector<complex> of size _N_pars
        inline std::vector<complex> complex_convert(const double * cpars)
        {
            std::vector<complex> result;
            // Move along the pars index when a parameter is not fixed
            int i = 0;
            for (auto par : _pars)
            {
                double mod = par._mod, arg = par._arg;
                if (!par._mod_fixed){ mod = cpars[i]; i++; };
                if (!par._arg_fixed){ arg = cpars[i]; i++; };
                result.push_back(mod*exp(I*arg));
            };

            if (i != _Nfree) warning("fitter::complex_convert", "Something went wrong in converting parameter vector.");
            return result;
        };

        
        // -----------------------------------------------------------------------
        // Methods to print out status to command line

        // Summary of data sets that have been recieved
        inline void data_info()
        {
            using std::cout; using std::left; using std::endl; using std::setw;
            
            if (_data.size() == 0){  warning("fitter::data_info", "No data found!"); return; };

            cout << left;
            divider();
            cout << "Fitting amplitude (\"" << _amplitude->name() << "\") to " << _N << " data points:" << endl;
            line();
            cout << setw(27) << "DATA SET"              << setw(25) << "TYPE     "               << setw(10) << "POINTS" << endl;
            cout << setw(27) << "---------------------" << setw(25) << "---------------------"   << setw(10) << "-------" << endl;
            for (auto data : _data)
            {
                cout << setw(27) << data._id  << setw(25)  << F::data_type(data._type)  << setw(10) << data._N << endl;  
            };
        };

        // Similar summary for parameters
        // Display alongside a vector of current parameter values
        // bool start is whether this is the starting guess vector or the 
        // best fit results
        inline void parameter_info()
        {
            using std::cout; using std::left; using std::endl; using std::setw;

            cout << std::setprecision(8);
            cout << left;

            line(); divider();
            // Print message at the beginning of the fit
            cout << "Fitting " + std::to_string(_Nfree) << " (of " << std::to_string(2*_pars.size()) << ") parameters" << endl;
            line();

            cout << left << setw(8) << "i"     << setw(15) << "PARAMETER"  << setw(26) << "START VALUE"       << endl;
            cout << left << setw(8) << "-----" << setw(15) << "----------" << setw(26) << "----------------------" << endl;
       
            for (auto par : _pars)
            {
                if (par._real) cout << left << setw(8) << par._i << setw(15) << par._label   << setw(26) << par._mod << setw(18) << par._mod_message << endl;
                else
                {
                    cout << left << setw(8) << par._i << setw(15) << par._label   << setw(26) << par.value() << endl;
                    cout << left << setw(8) << ""     << setw(15) << "  -> mod"   << setw(26) << par._mod << setw(18) << par._mod_message << endl;
                    cout << left << setw(8) << ""     << setw(15) << "  -> arg"   << setw(26) << par._arg << setw(18) << par._arg_message << endl;
                };
            };

            divider(); line();
        };

        // After a fit return a summary of fit results
        // At the end of a fit, print out a table sumarizing the fit results
        // if last_fit == true, we grab the results from the most recent fit in _minuit
        // else we print out the ones saved in _best_fit
        inline void print_results(bool last_fit = true)
        {
            using std::cout; using std::left; using std::endl; using std::setw;

            cout << std::setprecision(8);
            cout << left;

            double fcn               = _minuit->MinValue();
            double fcn_dof           = _minuit->MinValue() / double(dof());
            std::vector<double> pars = convert(_minuit->X());
            std::vector<double> errs = convert(_minuit->Errors());

            divider();
            std::cout << std::left << std::setw(5)  << "fcn = "      << std::setw(15) << fcn     << std::setw(5) << "";
            std::cout << std::left << std::setw(10) << "fcn/dof = "  << std::setw(15) << fcn_dof << "\n";

            line();

            cout << left << setw(8) << "i"     << setw(15) << "PARAMETER"  << setw(26) << "FIT VALUE"         << setw(18) << "ERROR"        << endl;
            cout << left << setw(8) << "-----" << setw(15) << "----------" << setw(26) << "----------------------" << setw(18) << "------------" << endl;

            int i = 0;
            for (auto par : _pars)
            {
                double mod = par._mod, arg = par._arg;
                std::string mod_err = par._mod_message, arg_err = par._arg_message;
                if (!par._mod_fixed){mod = pars[i]; mod_err = to_string(errs[i]); i++;};
                if (!par._arg_fixed){arg = pars[i]; arg_err = to_string(errs[i]); i++;};

                // Check arg if its between [-pi, pi]
                if (abs(arg) > PI) arg -= sign(arg)*2*PI;

                if (par._real) cout << left << setw(8) << par._i << setw(15) << par._label  << setw(26) << mod << setw(18) << mod_err << endl;
                else
                {
                    cout << left << setw(8) << par._i << setw(15) << par._label  << setw(26) << mod*exp(I*arg) << endl;
                    cout << left << setw(8) << ""     << setw(15) << "  -> mod"  << setw(26) << to_string(mod)  << setw(18) << mod_err << endl; 
                    cout << left << setw(8) << ""     << setw(15) << "  -> arg"  << setw(26) << to_string(arg)  << setw(18) << arg_err << endl; 
                };
            };
            divider(); line();
            
            // At the end update the amplitude parameters to include the fit results
            _amplitude->set_parameters(complex_convert(_minuit->X()));
        };
    };
}; // namespace iterateKT
#endif