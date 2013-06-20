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
## Configuration file for the c++/bin subdirectory
##
## author Come Raczy
##
################################################################################

include(${iSAAC_CXX_EXECUTABLE_CMAKE})

file (GLOB iSAAC_PROGRAM_SOURCE_LIST [a-zA-Z0-9]*.cpp)

##
## Generic rule for all the other programs
##
#set(BAM_LIBRARY "${SAMTOOLS_DIR}/libbam.a")
foreach(iSAAC_PROGRAM_SOURCE ${iSAAC_PROGRAM_SOURCE_LIST})
    get_filename_component(iSAAC_PROGRAM ${iSAAC_PROGRAM_SOURCE} NAME_WE)
    add_executable        (${iSAAC_PROGRAM} ${iSAAC_PROGRAM_SOURCE})
    target_link_libraries (${iSAAC_PROGRAM} ${iSAAC_AVAILABLE_LIBRARIES}
                           ${BAM_LIBRARY} ${Boost_LIBRARIES}
                           ${iSAAC_DEP_LIB} 
                           ${iSAAC_ADDITIONAL_LIB} )
    install(TARGETS ${iSAAC_PROGRAM} RUNTIME DESTINATION ${iSAAC_BINDIR})
endforeach(iSAAC_PROGRAM_SOURCE)
