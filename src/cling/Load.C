// Load script that links all the required libraries.
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@alumni.iu.edu
// -----------------------------------------------------------------------------

void Load()
{
    TString lib_ext   = gSystem->GetSoExt();

    //----------------------------------------------------------------------
    // Core library

    TString main_dir  = gSystem->Getenv("ITERATEDOKT");

    // Load the main library files
    TString lib  = main_dir + "/lib/libITERATEDOKT." + lib_ext;

    // Headers
    TString core = main_dir + "/src"; 
    TString phys = main_dir + "/models";
    TString data = main_dir + "/data";

    if (!gSystem->AccessPathName(lib.Data()))
    {
        Int_t lib_loaded = gSystem->Load(lib.Data());
        if (lib_loaded < 0) Fatal("Load()", "Library not loaded sucessfully!");

        gInterpreter->AddIncludePath( core.Data());
        gInterpreter->AddIncludePath( phys.Data());
        gInterpreter->AddIncludePath( data.Data());
    }
    else
    {
        Warning("Load()", "iteratedOKT library not found! Looked in: %s", lib.Data());
    }
}