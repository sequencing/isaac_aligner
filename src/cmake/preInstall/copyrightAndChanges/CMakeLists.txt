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
## Configuration file for the COPYRIGHT and Changes file installation
##
## author Roman Petrovski
##
################################################################################

include ("${iSAAC_GLOBALS_CMAKE}")

install(CODE "
    file(INSTALL DESTINATION \"${iSAAC_ORIG_DATADIR}\" TYPE FILE FILES \"${iSAAC_SOURCE_DIR}/Changes\" \"${CMAKE_SOURCE_DIR}/../COPYRIGHT\" PERMISSIONS ${iSAAC_LIBRARY_PERMISSIONS})
    ")
