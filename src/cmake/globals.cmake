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
## file globals.cmake
##
## CMake configuration file to identify the configuration of the system
##
## author Roman Petrovski
##
################################################################################

set(iSAAC_ORIG_DATADIR     "${iSAAC_PREFIX}${iSAAC_DATADIR}")
set(iSAAC_ORIG_BINDIR      "${iSAAC_EXEC_PREFIX}${iSAAC_BINDIR}")
set(iSAAC_ORIG_LIBDIR      "${iSAAC_EXEC_PREFIX}${iSAAC_LIBDIR}")
set(iSAAC_ORIG_LIBEXECDIR  "${iSAAC_EXEC_PREFIX}${iSAAC_LIBEXECDIR}")

#ensure that all file copying uses absolute paths
if    (NOT IS_ABSOLUTE "${iSAAC_PREFIX}")
    get_filename_component(iSAAC_ORIG_DATADIR     "${CMAKE_BINARY_DIR}/${iSAAC_ORIG_DATADIR}" ABSOLUTE)
endif (NOT IS_ABSOLUTE "${iSAAC_PREFIX}")

if    (NOT IS_ABSOLUTE "${iSAAC_EXEC_PREFIX}")
    get_filename_component(iSAAC_ORIG_BINDIR      "${CMAKE_BINARY_DIR}/${iSAAC_ORIG_BINDIR}" ABSOLUTE)
    get_filename_component(iSAAC_ORIG_LIBDIR      "${CMAKE_BINARY_DIR}/${iSAAC_ORIG_LIBDIR}" ABSOLUTE)
    get_filename_component(iSAAC_ORIG_LIBEXECDIR  "${CMAKE_BINARY_DIR}/${iSAAC_ORIG_LIBEXECDIR}" ABSOLUTE)
endif (NOT IS_ABSOLUTE "${iSAAC_EXEC_PREFIX}")

if    (iSAAC_MOVABLE_INSTALLATION)
# iSAAC_HOME must be set at this point to whatever would transform the component
# file path into the installation root path. for example for c++ bin module it is set to "../"
# Store relative component paths so that entire installation can be moved anywhere in the file system
        set(iSAAC_FULL_DATADIR      "${iSAAC_DATADIR}")
        set(iSAAC_FULL_BINDIR       "${iSAAC_BINDIR}")
        set(iSAAC_FULL_LIBDIR       "${iSAAC_LIBDIR}")
        set(iSAAC_FULL_LIBEXECDIR   "${iSAAC_LIBEXECDIR}")
else  (iSAAC_MOVABLE_INSTALLATION)
# Installation is not single-rooted. Use full paths for component resolution. Make sure iSAAC_HOME is empty.
        set(iSAAC_HOME       "")
        set(iSAAC_FULL_DATADIR      "${iSAAC_ORIG_DATADIR}")
        set(iSAAC_FULL_BINDIR       "${iSAAC_ORIG_BINDIR}")
        set(iSAAC_FULL_LIBDIR       "${iSAAC_ORIG_LIBDIR}")
        set(iSAAC_FULL_LIBEXECDIR   "${iSAAC_ORIG_LIBEXECDIR}")
endif (iSAAC_MOVABLE_INSTALLATION)

install(CODE "
    set(iSAAC_HOME      \"${iSAAC_HOME}\")
    set(iSAAC_FULL_DATADIR \"${iSAAC_FULL_DATADIR}\")
    set(iSAAC_FULL_BINDIR \"${iSAAC_FULL_BINDIR}\")
    set(iSAAC_FULL_LIBDIR \"${iSAAC_FULL_LIBDIR}\")
    set(iSAAC_FULL_LIBEXECDIR \"${iSAAC_FULL_LIBEXECDIR}\")
    ")

install(CODE "
    set(iSAAC_VERSION_FULL \"${iSAAC_VERSION_FULL}\")
    set(iSAAC_EXECUTABLE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
    set(iSAAC_LIBRARY_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
    ")

