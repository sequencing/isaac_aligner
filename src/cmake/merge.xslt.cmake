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
## Configuration file for the make subfolder
##
## author Roman Petrovski
##
################################################################################

# merge.xslt
add_custom_target(iSAAC_MERGE_XSLT
                   COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/${MERGE_XSLT}"
                   COMMAND tar -xzf "${MERGE_XSLT_REDIST_DIR}/${MERGE_XSLT}.tar.gz"
                   COMMAND mv "${CMAKE_CURRENT_BINARY_DIR}/${MERGE_XSLT}/${MERGE_XSLT}" "${CMAKE_CURRENT_BINARY_DIR}/${MERGE_XSLT}/MergeXmlDocuments.xsl"
                   COMMENT "Unpacking merge.xslt")

add_dependencies(iSAAC_OPT iSAAC_MERGE_XSLT)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${MERGE_XSLT}/MergeXmlDocuments.xsl" DESTINATION "${iSAAC_DATADIR}/xsl/common/")
