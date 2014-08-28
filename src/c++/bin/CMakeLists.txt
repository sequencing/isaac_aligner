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
## Configuration file for the c++/bin subdirectory
##
## author Come Raczy
##
################################################################################
set(iSAAC_HOME "../")
include(${iSAAC_CXX_EXECUTABLE_CMAKE})

file (GLOB iSAAC_PROGRAM_SOURCE_LIST [a-zA-Z0-9]*.cpp)


##
## Generic rule for all the other programs
##
foreach(iSAAC_PROGRAM_SOURCE ${iSAAC_PROGRAM_SOURCE_LIST})
    get_filename_component(iSAAC_FILENAME ${iSAAC_PROGRAM_SOURCE} NAME)
	configure_file(${iSAAC_PROGRAM_SOURCE} ${CMAKE_CURRENT_BINARY_DIR}/${iSAAC_FILENAME})
    get_filename_component(iSAAC_PROGRAM ${iSAAC_PROGRAM_SOURCE} NAME_WE)
    add_executable        (${iSAAC_PROGRAM} ${CMAKE_CURRENT_BINARY_DIR}/${iSAAC_FILENAME})
    target_link_libraries (${iSAAC_PROGRAM} ${iSAAC_AVAILABLE_LIBRARIES}
                           ${Boost_LIBRARIES} ${iSAAC_DEP_LIB}
                           ${iSAAC_ADDITIONAL_LIB} )
install(CODE "
    file(INSTALL DESTINATION \"${iSAAC_ORIG_BINDIR}\" TYPE FILE FILES \"${CMAKE_CURRENT_BINARY_DIR}/${iSAAC_PROGRAM}\" PERMISSIONS ${iSAAC_EXECUTABLE_PERMISSIONS})
  ")

endforeach(iSAAC_PROGRAM_SOURCE)
