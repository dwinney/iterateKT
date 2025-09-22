// Load script that links all the required libraries.
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@alumni.iu.edu
// -----------------------------------------------------------------------------

void Load()
{
    TString lib_ext   = gSystem->GetSoExt();
    
    // Fix for macOS - ROOT sometimes returns .so instead of .dylib
    #ifdef __APPLE__
    if (lib_ext == ".so") {
        lib_ext = ".dylib";
    }
    #endif

    //----------------------------------------------------------------------
    // Core library

    // Get the directory from the environment variable set by the executable
    TString main_dir = gSystem->Getenv("ITERATEKT_PATH");
    if (main_dir == "") {
        // Fallback: try to find the executable in the current working directory
        main_dir = gSystem->pwd();
    }

    // Load the main library files
    TString lib  = main_dir + "/lib/libITERATEKT." + lib_ext;

    // Headers
    TString core    = main_dir + "/src"; 
    TString physics = main_dir + "/physics";
    TString data    = main_dir + "/analysis";

    if (!gSystem->AccessPathName(lib.Data()))
    {
        Int_t lib_loaded = gSystem->Load(lib.Data());
        if (lib_loaded < 0) Fatal("Load()", "Library not loaded sucessfully!");

        gInterpreter->AddIncludePath( core.Data());
        gInterpreter->AddIncludePath( data.Data());
        gInterpreter->AddIncludePath( physics.Data());
        gInterpreter->AddIncludePath( main_dir.Data());
        
        // Add Boost include path
        gInterpreter->AddIncludePath("/opt/homebrew/include");
    }
    else
    {
        Warning("Load()", "iterateKT library not found! Looked in: %s", lib.Data());
    }
}