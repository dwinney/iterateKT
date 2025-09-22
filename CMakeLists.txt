cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.16)
cmake_policy(SET CMP0167 NEW)

project(iterateKT VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Validate C++ standard support
if(NOT CMAKE_CXX_STANDARD)
    message(FATAL_ERROR "C++17 standard is required but not supported by compiler ${CMAKE_CXX_COMPILER}")
endif()

# Detect operating system
if(APPLE)
    set(OS_NAME "macOS")
    set(OS_MACOS TRUE)
elseif(WIN32)
    set(OS_NAME "Windows")
    set(OS_WINDOWS TRUE)
elseif(UNIX)
    set(OS_NAME "Linux")
    set(OS_LINUX TRUE)
else()
    set(OS_NAME "Unknown")
endif()

message(STATUS "Detected OS: ${OS_NAME}")

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# Validate build type
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
    message(FATAL_ERROR "Invalid build type: ${CMAKE_BUILD_TYPE}. Valid options: Debug, Release, RelWithDebInfo, MinSizeRel")
endif()

# Set output directories (in build directory)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Platform-specific settings
if(OS_MACOS)
    # macOS-specific settings
    set(CMAKE_MACOSX_RPATH ON)
    set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    
    # Add -fPIC flag for macOS (avoid duplicates)
    if(NOT CMAKE_CXX_FLAGS MATCHES "-fPIC")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    endif()
    
    # Handle library naming for macOS
    set(CMAKE_SHARED_LIBRARY_PREFIX "lib")
    set(CMAKE_SHARED_LIBRARY_SUFFIX ".dylib")
    
elseif(OS_LINUX)
    # Linux-specific settings
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    
    # Add -fPIC flag for Linux (avoid duplicates)
    if(NOT CMAKE_CXX_FLAGS MATCHES "-fPIC")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    endif()
    
elseif(OS_WINDOWS)
    # Windows-specific settings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
    # Windows doesn't use RPATH, libraries must be in PATH or same directory
endif()

# Global RPATH settings for better cross-platform compatibility
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Set install directories (in project root for convenience)
set(CMAKE_INSTALL_BINDIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# Set the ITERATEKT path for the executable (point to project root)
set(ITERATEKT_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# Find required packages
find_package(PkgConfig QUIET)

##-----------------------------------------------------------------------
## CERN ROOT

# Try to find ROOT using multiple methods
find_package(ROOT QUIET)

if(NOT ROOT_FOUND)
    # Try using root-config
    find_program(ROOT_CONFIG_EXECUTABLE root-config)
    if(ROOT_CONFIG_EXECUTABLE)
        execute_process(
            COMMAND ${ROOT_CONFIG_EXECUTABLE} --prefix
            OUTPUT_VARIABLE ROOT_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(ROOT_PREFIX)
            list(APPEND CMAKE_MODULE_PATH "${ROOT_PREFIX}/share/root/cmake")
            find_package(ROOT REQUIRED COMPONENTS Core MathMore RIO Rint Hist Gpad Matrix)
        endif()
    endif()
endif()

if(ROOT_FOUND)
    message(STATUS "Found ROOT: ${ROOT_INCLUDE_DIRS}")
    
    # Check ROOT version
    if(ROOT_VERSION VERSION_LESS 6.0.0)
        message(FATAL_ERROR "ROOT 6.0.0+ required, found ${ROOT_VERSION}")
    endif()
    message(STATUS "✓ ROOT version ${ROOT_VERSION} (required: 6.0.0+)")
    
    include(${ROOT_USE_FILE})
    # ROOT_USE_FILE automatically configures CMAKE_CXX_FLAGS
else()
    message(FATAL_ERROR "ROOT not found! Please install ROOT first.\n"
        "Installation options:\n"
        "  - macOS: brew install root\n"
        "  - Linux: Use your package manager or build from source\n"
        "  - Windows: Download from https://root.cern/install/")
endif()

##-----------------------------------------------------------------------
## BOOST C++

set(BOOST_MIN_VERSION 1.68.0)
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED COMPONENTS system filesystem math)
if(Boost_FOUND)
    message(STATUS "✓ Boost version ${Boost_VERSION} (required: ${BOOST_MIN_VERSION}+)")
    message(STATUS "Found Boost: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost Libraries: ${Boost_LIBRARIES}")
else()
    message(FATAL_ERROR "Boost not found! Please install Boost first.\n"
        "Installation options:\n"
        "  - macOS: brew install boost\n"
        "  - Linux: sudo apt-get install libboost-all-dev (Ubuntu/Debian)\n"
        "  - Windows: Download from https://www.boost.org/users/download/")
endif()

##-----------------------------------------------------------------------
## Build and install

# Add source directory
add_subdirectory(src)

# Print configuration summary
message(STATUS "")
message(STATUS "Configuration Summary:")
message(STATUS "  OS: ${OS_NAME}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  C++ Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "  Install Prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  Install Bin Dir: ${CMAKE_INSTALL_BINDIR}")
message(STATUS "  Install Lib Dir: ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "")