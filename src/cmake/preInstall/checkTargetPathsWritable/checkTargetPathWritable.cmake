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
## file cppunit.cmake
##
## Configuration file for the cppunit subfolders
##
## author Roman Petrovski
##
################################################################################

foreach (iSAAC_DEST_DIR ${iSAAC_DEST_DIRS})
    message (STATUS "Testing access to ${iSAAC_DEST_DIR}...")
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "make_directory" "${iSAAC_DEST_DIR}/test" RESULT_VARIABLE TMP_RESULT )
    if (TMP_RESULT)
        message (STATUS "ERROR: Directory is not writeable: ${iSAAC_DEST_DIR}")
        message (STATUS "If you don't have administrator access to the "
                         "target installation location, please use --prefix "
                         "command-line option when configuring iSAAC. "
                         "Please use configure --help for all installer "
                         "command-line options details.")
        message (FATAL_ERROR "ERROR: iSAAC installation cannot continue")
    else (TMP_RESULT)
#        execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "remove_directory" "${iSAAC_DEST_DIR}/test" )
        message (STATUS "Directory is writeable: ${iSAAC_DEST_DIR}")
    endif (TMP_RESULT)
endforeach (iSAAC_DEST_DIR ${iSAAC_DEST_DIRS})
