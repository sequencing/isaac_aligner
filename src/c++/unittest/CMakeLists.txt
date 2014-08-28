################################################################################
##
## Isaac Genome Alignment Software
## Copyright (c) 2010-2014 Illumina, Inc.
## All rights reserved.
##
## This software is provided under the terms and conditions of the
## BSD 2-Clause License
##
## You should have received a copy of the BSD 2-Clause License
## along with this program. If not, see
## <https://github.com/sequencing/licenses/>.
##
################################################################################
##
## file CMakeLists.txt
##
## Configuration file for the c++/unittest subfolder
##
## author Come Raczy
##
################################################################################

include_directories(${CPPUNIT_INCLUDE_DIR})
add_library(isaac_cppunit cppunitTest.cpp RegistryName.cpp)
