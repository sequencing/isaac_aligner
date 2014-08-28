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
## Configuration file for the lib subfolder
##
## author Roman Petrovski
##
################################################################################

##
## List of iSAAC libraries
##
set (iSAAC_ALL_LIBRARIES
    bgzf 
    common 
    oligo 
    io
    rta
    flowcell 
    demultiplexing
    xml
    reference 
    bam 
    options 
    statistics 
    alignment 
    build 
    package
    reports 
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

