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

get_filename_component(iSAAC_CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message (STATUS "Adding the c++    program subdirectory: ${iSAAC_CURRENT_DIR_NAME}")
include_directories (${iSAAC_CXX_ALL_INCLUDES})
include_directories (${CMAKE_CURRENT_BINARY_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${iSAAC_CXX_CONFIG_H_DIR})

