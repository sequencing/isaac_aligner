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
## <https://github.com/downloads/sequencing/licenses/>.
##
## The distribution includes the code libraries listed below in the
## 'redist' sub-directory. These are distributed according to the
## licensing terms governing each library.
##
################################################################################
##
## file CMakeLists.txt
##
## Configuration file for the lib subfolder
##
## author Roman Petrovski
##
################################################################################

##
## List of iSAAC libraries
##
set (iSAAC_ALL_LIBRARIES 
    common 
    oligo 
    io
    basecalls
    flowcell 
    demultiplexing 
    reference 
    bam 
    options 
    statistics 
    alignment 
    build 
    reports 
    casava
    workflow)

##
## Build all the libraries for the project
## iSAAC_AVAILABLE_LIBRARIES is incrementally updated
##

set (iSAAC_AVAILABLE_LIBRARIES "")
set (iSAAC_ALL_LIBRARY_DIRS "")
foreach (iSAAC_LIB_DIR ${iSAAC_ALL_LIBRARIES})
    add_subdirectory(${iSAAC_LIB_DIR})
    set(iSAAC_AVAILABLE_LIBRARIES isaac_${iSAAC_LIB_DIR} ${iSAAC_AVAILABLE_LIBRARIES} )
endforeach (iSAAC_LIB_DIR)

set (iSAAC_AVAILABLE_LIBRARIES ${iSAAC_AVAILABLE_LIBRARIES} PARENT_SCOPE)

