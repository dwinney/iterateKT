# BUILD LIBRARY 
file(GLOB SRC "*.cpp" ) 

# Create the library target
add_library(ITERATEKT SHARED ${SRC} G__ITERATEKT.cxx )

# Set target-specific include directories instead of global ones
target_include_directories(ITERATEKT 
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${Boost_INCLUDE_DIRS}
)

ROOT_GENERATE_DICTIONARY(G__ITERATEKT
    LINKDEF ITERATEKTLinkDef.h
)

# Set library properties for better cross-platform compatibility
set_target_properties(ITERATEKT PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    OUTPUT_NAME "ITERATEKT"
)

# Link libraries with proper Boost components
target_link_libraries(ITERATEKT 
    ROOT::Core 
    ROOT::MathMore 
    ROOT::Matrix 
    ROOT::RIO 
    ROOT::Rint 
    ROOT::Gpad 
    ROOT::Hist
    ROOT::Minuit2
    Boost::system
    Boost::filesystem
    Boost::math
)

# Platform-specific post-processing
if(APPLE)
    # On macOS, create a .so symlink for ROOT compatibility
    add_custom_command(TARGET ITERATEKT POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
        $<TARGET_FILE_NAME:ITERATEKT>
        $<TARGET_FILE_DIR:ITERATEKT>/libITERATEKT.so
        COMMENT "Creating .so symlink for ROOT compatibility on macOS"
    )
    
    # Also create the symlink during install
    install(CODE "
        if(APPLE)
            execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                libITERATEKT.dylib
                ${CMAKE_INSTALL_LIBDIR}/libITERATEKT.so
            )
        endif()
    ")
endif()

# Install library
install(TARGETS ITERATEKT
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Install ROOT dictionary files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/../lib/libITERATEKT_rdict.pcm    
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/../lib/libITERATEKT.rootmap    
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

# BUILD EXECUTABLE
add_executable(iterateKT_exe "cling/iterateKT.cpp")

# Set executable properties
set_target_properties(iterateKT_exe PROPERTIES
    OUTPUT_NAME "iterateKT"
)

# Set target-specific include directories for executable
target_include_directories(iterateKT_exe 
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${Boost_INCLUDE_DIRS}
)

# Define the ITERATEKT path for the executable
target_compile_definitions(iterateKT_exe PRIVATE ITERATEKT_PATH="${ITERATEKT_PATH}")

# Link libraries
target_link_libraries(iterateKT_exe 
    ROOT::Core 
    ROOT::MathMore 
    ROOT::Matrix 
    ROOT::RIO 
    ROOT::Rint 
    ROOT::Minuit2
    ITERATEKT
)

# Install executable
install(TARGETS iterateKT_exe
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)