cmake_minimum_required (VERSION 3.30)
cmake_policy(VERSION 3.30)

project (iteratedKT)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-fPIC -O3") 
set(CMAKE_BUILD_TYPE "Release")

# Set output lib and bin locations
set(CMAKE_INSTALL_BINDIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_INSTALL_BINDIR})

# macOS is stupid, use these to get it to stop renaming install files
set(SKIP_BUILD_RPATH         TRUE)
set(CMAKE_SKIP_RPATH         TRUE)
set(CMAKE_SKIP_INSTALL_RPATH TRUE)
set(INSTALL_NAME_DIR           "")

##-----------------------------------------------------------------------
## CERN ROOT

# Make sure we can find ROOT 
execute_process(COMMAND root-config --prefix 
                COMMAND tr -d \\n
                RESULT_VARIABLE RC_RESULT 
                OUTPUT_VARIABLE ROOT_PREFIX)
list(APPEND CMAKE_MODULE_PATH "${ROOT_PREFIX}/share/root/cmake")

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Core MathMore RIO Rint Hist Gpad Matrix)
if (ROOT_FOUND)
    message("-- Found ROOT: ${ROOT_INCLUDE_DIRS}")
    #---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
    include(${ROOT_USE_FILE})
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
else()
    message(SEND_ERROR "-- Root not found!")
endif()

##-----------------------------------------------------------------------
## BOOST C++

set(BOOST_MIN_VERSION 1.68.0) # Need complex Gauss-Kronrod integration
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED)
if (Boost_FOUND)
    message("-- Found Boost: ${Boost_INCLUDE_DIRS}")
else()
    message(SEND_ERROR "-- Boost not found!")
endif()


##-----------------------------------------------------------------------
# Build core library 

add_subdirectory(src)