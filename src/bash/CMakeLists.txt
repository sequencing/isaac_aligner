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
## Configuration file for the bash subdirectory
##
## author Roman Petrovski
##
################################################################################
set(iSAAC_HOME "\\\$(dirname \\\$(readlink -e \\\$0))/../")

add_subdirectory (bin)
add_subdirectory (libexec)

