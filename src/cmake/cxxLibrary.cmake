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
## file cxxLibrary.cmake
##
## CMake configuration file for all the c++ libraries
##
## author Come Raczy
##
################################################################################
include_directories (${iSAAC_CXX_ALL_INCLUDES})
include_directories (${CMAKE_CURRENT_BINARY_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${iSAAC_CXX_CONFIG_H_DIR})

get_filename_component(iSAAC_CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message (STATUS "Adding the c++    library subdirectory: ${iSAAC_CURRENT_DIR_NAME}")

##
## Some generators (VS) require all targets to be unique across the project.
## Therefore, a unique prefix is needed to create the target names which are
## shared across libraries
##

string(REGEX REPLACE ${CMAKE_SOURCE_DIR}/c[+][+]/ "" TMP1 ${CMAKE_CURRENT_SOURCE_DIR}/)
string(REGEX REPLACE "/" "_" iSAAC_UNIQUE_PREFIX ${TMP1})

##
## build the library
##
file(GLOB_RECURSE iSAAC_LIBRARY_SOURCES_WITH_CPPUNIT *.cpp *.c)

foreach (SOURCE_FILE ${iSAAC_LIBRARY_SOURCES_WITH_CPPUNIT})
    string(REGEX MATCH "cppunit" CPPUNIT_MATCH ${SOURCE_FILE} )
    if (NOT CPPUNIT_MATCH)
        set(iSAAC_LIBRARY_SOURCES ${iSAAC_LIBRARY_SOURCES} ${SOURCE_FILE})
    endif (NOT CPPUNIT_MATCH)
endforeach (SOURCE_FILE)

foreach (SOURCE_FILE ${iSAAC_LIBRARY_SOURCES})
    get_filename_component(SOURCE_NAME ${SOURCE_FILE} NAME_WE)
    if (${SOURCE_NAME}_COMPILE_FLAGS)
        set_source_files_properties(${SOURCE_FILE} PROPERTIES COMPILE_FLAGS ${${SOURCE_NAME}_COMPILE_FLAGS})
    endif (${SOURCE_NAME}_COMPILE_FLAGS)
endforeach (SOURCE_FILE)

#include_directories (${iSAAC_COMMON_INCLUDE} )
add_library         (isaac_${iSAAC_LIB_DIR} STATIC ${iSAAC_LIBRARY_SOURCES})
add_dependencies(isaac_${iSAAC_LIB_DIR} iSAAC_OPT)

##
## build the unit tests if any (this should be mandatory really)
##

if (HAVE_CPPUNIT AND iSAAC_UNIT_TESTS)
    find_path(${CMAKE_CURRENT_SOURCE_DIR}_CPPUNIT_DIR cppunit PATHS ${CMAKE_CURRENT_SOURCE_DIR} NO_DEFAULT_PATH)
    if (${CMAKE_CURRENT_SOURCE_DIR}_CPPUNIT_DIR)
        message (STATUS "Adding the cppunit subdirectory for ${iSAAC_LIB_DIR}")
        add_subdirectory (cppunit)
    endif (${CMAKE_CURRENT_SOURCE_DIR}_CPPUNIT_DIR)
endif(HAVE_CPPUNIT AND iSAAC_UNIT_TESTS)

