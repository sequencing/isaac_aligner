################################################################################
##
## Isaac Genome Alignment Software
## Copyright (c) 2010-2012 Illumina, Inc.
##
## This software is provided under the terms and conditions of the
## Illumina Open Source Software License 1.
##
## You should have received a copy of the Illumina Open Source
## Software License 1 along with this program. If not, see
## <https://github.com/sequencing/licenses/>.
##
## The distribution includes the code libraries listed below in the
## 'redist' sub-directory. These are distributed according to the
## licensing terms governing each library.
##
################################################################################
##
## file CMakeLists.txt
##
## Configuration file for the c++ subfolder
##
## author Come Raczy
##
################################################################################
include (${iSAAC_GLOBALS_CMAKE})

set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set (iSAAC_CXX_EXECUTABLE_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxExecutable.cmake")
set (iSAAC_CXX_LIBRARY_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxLibrary.cmake")
set (iSAAC_CXX_CONFIGURE_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxConfigure.cmake")

set(iSAAC_CXX_CONFIG_H_DIR ${CMAKE_CURRENT_BINARY_DIR}/common)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
include ("${iSAAC_CXX_CONFIGURE_CMAKE}")

if (HAVE_CPPUNIT AND iSAAC_UNIT_TESTS)
    set (iSAAC_CPPUNIT_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cppunit.cmake")
    add_subdirectory (unittest)
endif (HAVE_CPPUNIT AND iSAAC_UNIT_TESTS)


##
## The include directories
##
set (iSAAC_CXX_ALL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_BINARY_DIR}")

##
## Build all the libraries for the project
##
add_subdirectory (lib)

##
## build all the applications for the project
##

add_subdirectory (bin)
add_subdirectory (libexec)

##
## build all the internal applications for the project
##

#add_subdirectory (libexec)

##
## build the documentation when available
##
include  (FindDoxygen)
message (STATUS "Doxygen: ${DOXYGEN_EXECUTABLE}. Dot: ${DOXYGEN_DOT_EXECUTABLE}.")
if (DOXYGEN_FOUND)
    set (DOXYFILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    message (STATUS "Creating Doxygen config file: ${DOXYFILE}")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${DOXYFILE} @ONLY IMMEDIATE)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOXYFILE})
endif (DOXYGEN_FOUND)
