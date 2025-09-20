// Buildable executable that loads Load.C and runs a script
//
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// -----------------------------------------------------------------------------

#include <iostream>
#include <TRint.h>
#include <TEnv.h>
#include <TString.h>
#include <TSystem.h>
#include <vector>

int main(int argc, char **argv)
{
    TString macroName;
    for(Int_t i = 0; i < argc; i++)
    {
        TString opt = argv[i];
        if((opt.Contains(".C")))   macroName = opt;
        if((opt.Contains(".cpp"))) macroName = opt;
    }

    // All of this is just to add -l and -b to argv so that 
    // ROOT opens without a welcome message and in batch mode...
    std::vector<char*> new_argv(argv, argv + argc);
    new_argv.push_back((char *) "-l");
    new_argv.push_back((char *) "-b");
    new_argv.push_back(NULL);
    argv = new_argv.data();
    argc += 2;    

    TRint * app = new TRint( "iterateKT", &argc, argv);
    
    // Use the predefined path from CMake
    TString env = ITERATEKT_PATH;
    
    // Set a ROOT variable that Load.C can use
    app->ProcessLine(Form("gSystem->Setenv(\"ITERATEKT_PATH\", \"%s\");", env.Data()));
    
    app->ProcessLine(Form(".x %s/src/cling/Load.C", env.Data()));
    app->ProcessLine(Form(".x %s", macroName.Data()));
    app->Terminate(0);

    delete app;
    return 0;
};