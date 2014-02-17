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
## file cxxExecutable.cmake
##
## CMake configuration file for all the c++ executables
##
## author Roman Petrovski
##
################################################################################

include (${iSAAC_GLOBALS_CMAKE})

# Support for static linking. Notice this is done here and not in cxxConfigure to 
# allow dynamic linking for cppunit tests as some platforms lack static libcppunit.
# Note that this implies that all libraries must be found with the
# exact file name (libXXX.a or libXXX.so)
if    (iSAAC_FORCE_STATIC_LINK)
    message(STATUS "All libraries will be statically linked")
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-static")
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-static")
    # ensure that even if cmake decides to allow for dynamic libs resolution,
    # this gets overriden into static...
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS ${CMAKE_EXE_LINK_STATIC_CXX_FLAGS})
    set(iSAAC_LIBRARY_PREFIX ${CMAKE_STATIC_LIBRARY_PREFIX})
    set(iSAAC_LIBRARY_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
else  (iSAAC_FORCE_STATIC_LINK)
    set(iSAAC_LIBRARY_PREFIX "")
    set(iSAAC_LIBRARY_SUFFIX "")
endif (iSAAC_FORCE_STATIC_LINK)

get_filename_component(iSAAC_CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message (STATUS "Adding the c++    program subdirectory: ${iSAAC_CURRENT_DIR_NAME}")
include_directories (${iSAAC_CXX_ALL_INCLUDES})
include_directories (${CMAKE_CURRENT_BINARY_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${iSAAC_CXX_CONFIG_H_DIR})

