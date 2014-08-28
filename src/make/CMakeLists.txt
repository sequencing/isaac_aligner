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
## Configuration file for the make subfolder
##
## author Roman Petrovski
##
################################################################################
set(iSAAC_HOME "\$(dir \$(lastword \$(MAKEFILE_LIST)))../../../../")
include ("${iSAAC_GLOBALS_CMAKE}")

install(CODE "
    include (\"${iSAAC_MACROS_CMAKE}\")
    configure_files_recursively (\"${CMAKE_CURRENT_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}\" \"*.in\")
    install_files_recursively (\"${CMAKE_CURRENT_BINARY_DIR}\" \"${iSAAC_ORIG_DATADIR}/makefiles\" \"*.in\" \"\${iSAAC_LIBRARY_PERMISSIONS}\")
    configure_files_recursively (\"${CMAKE_CURRENT_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}\" \"*.mk\")
    install_files_recursively (\"${CMAKE_CURRENT_BINARY_DIR}\" \"${iSAAC_ORIG_DATADIR}/makefiles\" \"*.mk\" \"\${iSAAC_LIBRARY_PERMISSIONS}\")
    ")
