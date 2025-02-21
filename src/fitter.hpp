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
        : _fixed(true), 
          _mod(1.), _arg(0.),
          _label(default_label(0)), _i(-1)
        {};
        
        // Constructor for amplitude variables
        parameter(int i) : _i(i), _label(default_label(i)){};

        uint         _i;
        std::string _label;
        std::string _message;
        bool   _fixed         = false;
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
        // Methods to add data to be fit against

        inline void add_data(data_set data){ _N += data._N; _data.push_back(data); };
        inline void add_data(std::vector<data_set> data){ for (auto datum : data) add_data(datum); };
        inline void clear_data(){ _data.clear(); _N = 0; };

        // -----------------------------------------------------------------------
        // Set limits, labels, and fix parameters

        // Reset labels, limits and options on all parameters
        inline void reset_parameters()
        {
            std::vector<std::string> labels; // Default label names
            _pars.clear(); _Nfree = _amplitude->N_pars();
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

        inline void fix_parameter(parameter& par, complex val)
        {
            // If parameter is already fixed, just update the fixed val
            // otherwise flip the fixed flag and update the number of free pars
            if (!par._fixed) _Nfree--;
            par._fixed   = true;
            par._mod     = std::abs(val);
            par._arg     = std::arg(val);
            par._message = "[FIXED]";
        };

        inline void fix_parameter(std::string label, complex val)
        {
            int index = find_parameter(label);
            if (index < 0) return;
            return fix_parameter(_pars[index], val);
        };

        inline void free_parameter(parameter& par)
        {
            // if not fixed, this does nothing
            if (!par._fixed) return;
            par._fixed   = false;
            par._message = "";
            _Nfree++;
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
            if (starting_guess.size() != _Nfree) 
            {
                warning("fitter::do_fit", "Starting guess not the correct size! Expected " + std::to_string(_Nfree) + " parameters!");
                return;
            };

            set_up(starting_guess);

            if (show_data) { line(); data_info(); };
            parameter_info(starting_guess);

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
        inline void set_max_calls(int n){ print(n); _max_calls = n; };
        
        // Message level for minuit (0-4)
        inline void set_print_level(int n){ _print_level = n; };

        // Change tolerance
        inline void set_tolerance(double tol){ _tolerance = tol; };

        // private:

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
            int i = 0;
            for (auto &par : _pars)
            {   
                if (par._fixed) continue;
                
                // Set up the starting guess
                par._mod = std::abs(starting_guess[i]);
                par._arg = std::arg(starting_guess[i]);
                _minuit->SetVariable(i,   "|" + par._label + "|",    std::abs(starting_guess[i]), par._step);
                _minuit->SetVariableLowerLimit(i, 0.);    // mod is positive definite
                
                _minuit->SetVariable(i+1, "arg(" + par._label + ")", std::arg(starting_guess[i]), par._step);
                _minuit->SetVariableLimits(i+1, -PI, PI); // angle from [-pi,pi]

                i+=2; // move index up
            };
        
            _wfcn = ROOT::Math::Functor(this, &fitter::fit_fcn, 2*_Nfree);
            _minuit->SetFunction(_wfcn);
        };

                
        // -----------------------------------------------------------------------
        // Calcualtions of chi-squared      
      
        // This is the actual function that gets called by minuit
        inline double fit_fcn(const double *cpars)
        { 
            // First convert the C string to a C++ vector
            std::vector<complex> pars = complex_convert(cpars);
            // Pass parameters to the amplitude
            _amplitude->set_parameters(pars);

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
                if (par._fixed) result.push_back(par.value());
                else { result.push_back(cpars[i]*exp(I*cpars[i+1])); i+=2; };
            };

            if (i != 2*_Nfree) warning("fitter::complex_convert", "Something went wrong in converting parameter vector.");
            return result;
        };

        
        // -----------------------------------------------------------------------
        // Methods to print out status to command line

        // Summary of data sets that have been recieved
        inline void data_info()
        {
            using std::cout; using std::left; using std::endl; using std::setw;
            
            if (_data.size() == 0)
            {
                warning("fitter::data_info", "No data found!"); 
                return;
            };

            cout << left;
            divider();
            cout << "Fitting amplitude (\"" << _amplitude->name() << "\") to " << _N << " data points:" << endl;
            line();
            cout << setw(25) << "DATA SET"         << setw(30) << "TYPE     "      << setw(10) << "POINTS" << endl;
            cout << setw(25) << "----------------" << setw(30) << "--------------" << setw(10) << "-------" << endl;
            for (auto data : _data)
            {
                cout << setw(25) << data._id  << setw(30)  << F::data_type(data._type)  << setw(10) << data._N << endl;  
            };
        };

        // Similar summary for parameters
        // Display alongside a vector of current parameter values
        // bool start is whether this is the starting guess vector or the 
        // best fit results
        inline void parameter_info(std::vector<complex> starting_guess)
        {
            using std::cout; using std::left; using std::endl; using std::setw;

            cout << std::setprecision(8);
            cout << left;

            line(); divider();
            // Print message at the beginning of the fit
            cout << "Fitting " + std::to_string(_Nfree) << " (of " << std::to_string(_pars.size()) << ") parameters" << endl;
            line();

            cout << left << setw(10) << "i"     << setw(17) << "PARAMETER"  << setw(20) << "START VALUE"  << endl;
            cout << left << setw(10) << "-----" << setw(17) << "----------" << setw(20) << "------------" << endl;
       
            for (auto par : _pars)
            {
                cout << left << setw(10) << par._i << setw(17) << par._label  << setw(20) << endl;
                cout << left << setw(10) << ""     << setw(17) << "  -> mod"   << setw(20) << par._mod << setw(20) << par._message << endl;
                cout << left << setw(10) << ""     << setw(17) << "  -> arg"   << setw(20) << par._arg << setw(20) << par._message << endl;
            };

            line(); divider(); line();
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

            int dof                  = _N - _minuit->NFree();
            double fcn               = _minuit->MinValue();
            double fcn_dof           = _minuit->MinValue() / double(dof);
            std::vector<double> pars = convert(_minuit->X());
            std::vector<double> errs = convert(_minuit->Errors());
            
            divider();
            std::cout << std::left << std::setw(5)  << "fcn = "      << std::setw(15) << fcn     << std::setw(5) << "";
            std::cout << std::left << std::setw(10) << "fcn/dof = "  << std::setw(15) << fcn_dof << "\n";

            line();

            cout << left << setw(10) << "i"     << setw(16) << "PARAMETER"  << setw(18) << "FIT VALUE"    << setw(18) << "ERROR"        << endl;
            cout << left << setw(10) << "-----" << setw(16) << "----------" << setw(18) << "------------" << setw(18) << "------------" << endl;

            int i = 0;
            for (auto par : _pars)
            {
                if (par._fixed)
                {
                    cout << left << setw(10) << par._i << setw(16) << par._label  << setw(18) << endl;
                    cout << left << setw(10) << ""     << setw(16) << "  -> mod"   << setw(18) << par._mod << setw(18) << par._message << endl;
                    cout << left << setw(10) << ""     << setw(16) << "  -> arg"   << setw(18) << par._arg << setw(18) << par._message << endl;
                    continue;
                }

                cout << left << setw(10) << par._i << setw(16) << par._label  << setw(18) << endl;
                cout << left << setw(10) << ""     << setw(16) << "  -> mod"   << setw(18) << to_string(pars[i])   << setw(18) << to_string(errs[i])   << endl;
                cout << left << setw(10) << ""     << setw(16) << "  -> arg"   << setw(18) << to_string(pars[i+1]) << setw(18) << to_string(errs[i+1]) << endl;
                i+=2;
            };
            line(); divider(); line();
            
            // At the end update the amplitude parameters to include the fit results
            _amplitude->set_parameters(complex_convert(_minuit->X()));
        };
    };
}; // namespace iterateKT
#endif