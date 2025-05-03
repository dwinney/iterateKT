# BUILD LIBRARY 
file(GLOB SRC "*.cpp" ) 
include_directories($ENV{ITERATEKT}/src)

ROOT_GENERATE_DICTIONARY(G__ITERATEKT
    LINKDEF ITERATEKTLinkDef.h
)
add_library(ITERATEKT SHARED ${SRC} G__ITERATEKT.cxx )
target_link_libraries( ITERATEKT ROOT::Core ROOT::MathMore ROOT::Matrix ROOT::RIO ROOT::Rint ROOT::Gpad ROOT::Hist)
target_link_libraries( ITERATEKT Boost::boost)

install( TARGETS ITERATEKT
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )
install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libITERATEKT_rdict.pcm    
         DESTINATION ${CMAKE_INSTALL_LIBDIR})

# BUILD EXECUTABLE
add_executable( iterateKT "cling/iterateKT.cpp" )
target_link_libraries( iterateKT ROOT::Core ROOT::MathMore ROOT::Matrix ROOT::RIO ROOT::Rint )
target_link_libraries( iterateKT ITERATEKT )