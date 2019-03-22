#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
set(DK2NU $ENV{DK2NU})
set(ROOTSYS $ENV{ROOTSYS})
set(FHICLCPP_INC $ENV{FHICLCPP_INC})
set(CETLIB_EXCEPT_INC $ENV{CETLIB_EXCEPT_INC})
set(CETLIB_INC $ENV{CETLIB_INC})
set(CLHEP_INC $ENV{CLHEP_INC})

include_directories(${PROJECT_SOURCE_DIR}
		    ${FHICLCPP_INC}
		    ${CETLIB_EXCEPT_INC}
		    ${CETLIB_INC}
		    ${CLHEP_INC}
		    ${DK2NU}/include
		    ${ROOTSYS}/include )	    

file(GLOB sources  ${PROJECT_SOURCE_DIR}/Base/*.cxx)

add_library(BNBSysBase SHARED ${sources})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS BNBSysBase  DESTINATION ${PROJECT_SOURCE_DIR}/lib)

