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
## file configure
##
## Configuration file for the share subdirectory
##
## author Roman Petrovski
##
################################################################################

#include ("${iSAAC_MACROS_CMAKE}")

#get_filename_component(iSAAC_CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message (STATUS "Verifying target directories access")

include ("${iSAAC_GLOBALS_CMAKE}")
install(
    CODE "set(iSAAC_DEST_DIRS \"\${iSAAC_DEST_ETCDIR}\" \"\${iSAAC_DEST_DATADIR}\" 
                               \"\${iSAAC_DEST_BINDIR}\" \"\${iSAAC_DEST_LIBDIR}\" 
                               \"\${iSAAC_DEST_LIBEXECDIR}\")"
    SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/checkTargetPathWritable.cmake"
)
