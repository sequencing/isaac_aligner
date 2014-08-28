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
## Configuration file for the markdown subfolder
##
## author Roman Petrovski
##
################################################################################
include ("${iSAAC_GLOBALS_CMAKE}")

##
## build the documentation when available
##
include(FindMarkdown)
include(UseMarkdown)
message (STATUS "Markdown: ${MARKDOWN_EXECUTABLE}.")
if (MARKDOWN_FOUND)
    set(MARKDOWN_PARAMS -x tables -x toc)
    add_markdown_target(manual ${CMAKE_CURRENT_BINARY_DIR} "manual.md" RENAME "manual.html")
install(CODE "
    file(INSTALL DESTINATION \"${iSAAC_ORIG_DATADIR}\" TYPE FILE FILES \"${CMAKE_CURRENT_BINARY_DIR}/manual.html\" PERMISSIONS ${iSAAC_LIBRARY_PERMISSIONS})
  ")

#    install_markdown_target(manual DESTINATION \"${iSAAC_ORIG_DATADIR}/html\")

endif (MARKDOWN_FOUND)


