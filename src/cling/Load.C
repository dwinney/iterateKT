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

    TString main_dir  = gSystem->Getenv("ITERATEKT");

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
    }
    else
    {
        Warning("Load()", "iterateKT library not found! Looked in: %s", lib.Data());
    }

    // Add JSON library
    TString json = gSystem->Getenv("JSON");
    gInterpreter->AddIncludePath(json.Data());

    // Plotting library
    TString plot_dir = gSystem->Getenv("JPACUTILS");
    TString plot_lib = plot_dir+"/lib/libJPACUTILS."+lib_ext;
    TString plot_inc = plot_dir;

    if (!gSystem->AccessPathName(plot_lib.Data()))
    {
        gSystem->Load(plot_lib.Data());
        gInterpreter->AddIncludePath(plot_inc.Data());
    };
}