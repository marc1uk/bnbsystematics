#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
set(DK2NU $ENV{DK2NU})
set(ROOTSYS $ENV{ROOTSYS})
set(FHICLCPP_INC $ENV{FHICLCPP_INC})
set(CETLIB_EXCEPT_INC $ENV{CETLIB_EXCEPT_INC})
set(CETLIB_INC $ENV{CETLIB_INC})
set(CLHEP_INC $ENV{CLHEP_INC})
set(BOOSTLIB "-L$ENV{BOOST_LIB} -lboost_program_options -lboost_filesystem")

include_directories(${PROJECT_SOURCE_DIR}
		    ${FHICLCPP_INC}
		    ${CETLIB_EXCEPT_INC}
		    ${CETLIB_INC}
		    ${CLHEP_INC}
		    ${DK2NU}/include
		    ${ROOTSYS}/include)
#include BOOST as system (-isystem instead of -I) to suppress compiler warnings
include_directories(SYSTEM $ENV{BOOST_INC})


add_executable(rwgh rwgh.cxx)
find_library(fhiclcpp_LIBRARY NAMES fhiclcpp PATHS $ENV{FHICLCPP_LIB})
target_link_libraries(rwgh BNBSysBase BNBSysCalc ${ROOT_LIBRARIES} ${CLHEP_LIBRARIES} ${fhiclcpp_LIBRARY} ${BOOSTLIB} ${DK2NU}/lib/libdk2nuTree.so)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS rwgh  DESTINATION ${PROJECT_SOURCE_DIR}/bin)

