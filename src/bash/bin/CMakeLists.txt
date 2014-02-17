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
## Configuration file for the bin subfolder
##
## author Roman Petrovski
##
################################################################################
set(iSAAC_HOME "\\\$(dirname \\\$(readlink -e \\\$0))/../")
include ("${iSAAC_GLOBALS_CMAKE}")

install(CODE "
    include (\"${iSAAC_MACROS_CMAKE}\")
    configure_files_recursively (\"${CMAKE_CURRENT_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}\" \"isaac-*\")
    install_files_recursively (\"${CMAKE_CURRENT_BINARY_DIR}\" \"${iSAAC_ORIG_BINDIR}\" \"isaac-*\" \"\${iSAAC_EXECUTABLE_PERMISSIONS}\")
    ")

