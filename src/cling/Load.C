// Load script that links all the required libraries.
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@alumni.iu.edu
// -----------------------------------------------------------------------------

void Load()
{
    // Determine possible library file extensions
    TString so_ext   = gSystem->GetSoExt();
    TString main_dir = gSystem->Getenv("ITERATEKT");

    // Candidate library paths (try detected ext first, then the common alternative)
    std::vector<TString> candidates;
    candidates.push_back(main_dir + "/lib/libITERATEKT." + so_ext);
#ifdef __APPLE__
    if (so_ext != "dylib") candidates.push_back(main_dir + "/lib/libITERATEKT.dylib");
#else
    if (so_ext != "so")    candidates.push_back(main_dir + "/lib/libITERATEKT.so");
#endif

    // Headers (always add include paths so interpreted headers work even if lib load fails)
    TString core    = main_dir + "/src"; 
    TString physics = main_dir + "/physics";
    TString data    = main_dir + "/analysis";
    gInterpreter->AddIncludePath( core.Data());
    gInterpreter->AddIncludePath( data.Data());
    gInterpreter->AddIncludePath( physics.Data());
    gInterpreter->AddIncludePath( main_dir.Data());

    // Add Boost include path so interpreted headers can find <boost/...>
    TString boost_inc = gSystem->Getenv("BOOST_INCLUDEDIR");
    if (boost_inc.Length() == 0)
    {
        TString boost_root = gSystem->Getenv("BOOST_ROOT");
        if (boost_root.Length() > 0)
        {
            boost_inc = boost_root + "/include";
        }
#ifdef __APPLE__
        if (boost_inc.Length() == 0)
        {
            // Try Homebrew prefix
            TString brew_prefix = gSystem->GetFromPipe("brew --prefix boost 2>/dev/null");
            brew_prefix = brew_prefix.Strip(TString::kBoth, '\n');
            if (brew_prefix.Length() > 0) boost_inc = brew_prefix + "/include";
        }
#endif
    }
    if (boost_inc.Length() > 0) gInterpreter->AddIncludePath(boost_inc.Data());

    // Try to load the first library that exists
    bool loaded = false;
    for (auto &path : candidates)
    {
        if (!gSystem->AccessPathName(path.Data()))
        {
            Int_t lib_loaded = gSystem->Load(path.Data());
            if (lib_loaded < 0) Fatal("Load()", "Library not loaded sucessfully! Tried: %s", path.Data());
            loaded = true;
            break;
        }
    }
    if (!loaded)
    {
        Warning("Load()", "iterateKT library not found! Looked in: %s and fallbacks.", candidates.front().Data());
    }
}