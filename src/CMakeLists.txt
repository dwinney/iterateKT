# BUILD LIBRARY 
file(GLOB SRC "*.cpp" ) 
# Exclude the executable source from the library sources
list(REMOVE_ITEM SRC "${CMAKE_CURRENT_SOURCE_DIR}/cling/iterateKT.cpp")
include_directories($ENV{ITERATEKT}/src)

ROOT_GENERATE_DICTIONARY(G__ITERATEKT
    LINKDEF ITERATEKTLinkDef.h
)
add_library(ITERATEKT SHARED ${SRC} G__ITERATEKT.cxx )
target_link_libraries( ITERATEKT ROOT::Core ROOT::MathMore ROOT::Matrix ROOT::RIO ROOT::Rint ROOT::Gpad ROOT::Hist)
target_link_libraries( ITERATEKT Boost::boost)

set_target_properties(ITERATEKT PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_RPATH}"
    BUILD_RPATH   "${CMAKE_BUILD_RPATH}"
)

install( TARGETS ITERATEKT
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )
install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libITERATEKT_rdict.pcm    
         DESTINATION ${CMAKE_INSTALL_LIBDIR})

# BUILD EXECUTABLE (use different target name to avoid case-collision on macOS)
add_executable( iterateKT_app "cling/iterateKT.cpp" )
target_link_libraries( iterateKT_app ROOT::Core ROOT::MathMore ROOT::Matrix ROOT::RIO ROOT::Rint )
target_link_libraries( iterateKT_app ITERATEKT )

set_target_properties(iterateKT_app PROPERTIES
    OUTPUT_NAME "iterateKT"
    INSTALL_RPATH "${CMAKE_INSTALL_RPATH}"
    BUILD_RPATH   "${CMAKE_BUILD_RPATH}"
)