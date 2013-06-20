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
## file globals.cmake
##
## CMake configuration file to identify the configuration of the system
##
## author Roman Petrovski
##
################################################################################

set(iSAAC_ORIG_ETCDIR      "${CMAKE_INSTALL_PREFIX}/${iSAAC_ETCDIR}")
set(iSAAC_ORIG_DATADIR     "${CMAKE_INSTALL_PREFIX}/${iSAAC_DATADIR}")
set(iSAAC_ORIG_BINDIR      "${CMAKE_INSTALL_PREFIX}/${iSAAC_BINDIR}")
set(iSAAC_ORIG_LIBDIR      "${CMAKE_INSTALL_PREFIX}/${iSAAC_LIBDIR}")
set(iSAAC_ORIG_LIBEXECDIR  "${CMAKE_INSTALL_PREFIX}/${iSAAC_LIBEXECDIR}")

get_filename_component(iSAAC_FULL_ETCDIR       "${iSAAC_ORIG_ETCDIR}" ABSOLUTE)
get_filename_component(iSAAC_FULL_DATADIR      "${iSAAC_ORIG_DATADIR}" ABSOLUTE)
get_filename_component(iSAAC_FULL_BINDIR       "${iSAAC_ORIG_BINDIR}" ABSOLUTE)
get_filename_component(iSAAC_FULL_LIBDIR       "${iSAAC_ORIG_LIBDIR}" ABSOLUTE)
get_filename_component(iSAAC_FULL_LIBEXECDIR   "${iSAAC_ORIG_LIBEXECDIR}" ABSOLUTE)

set(iSAAC_PARTIAL_ETCDIR "${iSAAC_ETCDIR}")
set(iSAAC_PARTIAL_DATADIR "${iSAAC_DATADIR}")
set(iSAAC_PARTIAL_BINDIR "${iSAAC_BINDIR}")
set(iSAAC_PARTIAL_LIBDIR "${iSAAC_LIBDIR}")
set(iSAAC_PARTIAL_LIBEXECDIR "${iSAAC_LIBEXECDIR}")

install(CODE "

    # _DEST_ variables always point to location where files are copied
    get_filename_component(iSAAC_DEST_ETCDIR       \"\$ENV{DESTDIR}${iSAAC_ORIG_ETCDIR}\" ABSOLUTE)
    get_filename_component(iSAAC_DEST_DATADIR      \"\$ENV{DESTDIR}${iSAAC_ORIG_DATADIR}\" ABSOLUTE)
    get_filename_component(iSAAC_DEST_BINDIR       \"\$ENV{DESTDIR}${iSAAC_ORIG_BINDIR}\" ABSOLUTE)
    get_filename_component(iSAAC_DEST_LIBDIR       \"\$ENV{DESTDIR}${iSAAC_ORIG_LIBDIR}\" ABSOLUTE)
    get_filename_component(iSAAC_DEST_LIBEXECDIR   \"\$ENV{DESTDIR}${iSAAC_ORIG_LIBEXECDIR}\" ABSOLUTE)

    set(iSAAC_FULL_ETCDIR \"${iSAAC_FULL_ETCDIR}\")
    set(iSAAC_FULL_DATADIR \"${iSAAC_FULL_DATADIR}\")
    set(iSAAC_FULL_BINDIR \"${iSAAC_FULL_BINDIR}\")
    set(iSAAC_FULL_LIBDIR \"${iSAAC_FULL_LIBDIR}\")
    set(iSAAC_FULL_LIBEXECDIR \"${iSAAC_FULL_LIBEXECDIR}\")

    set(iSAAC_PARTIAL_ETCDIR \"${iSAAC_PARTIAL_ETCDIR}\")
    set(iSAAC_PARTIAL_DATADIR \"${iSAAC_PARTIAL_DATADIR}\")
    set(iSAAC_PARTIAL_BINDIR \"${iSAAC_PARTIAL_BINDIR}\")
    set(iSAAC_PARTIAL_LIBDIR \"${iSAAC_PARTIAL_LIBDIR}\")
    set(iSAAC_PARTIAL_LIBEXECDIR \"${iSAAC_PARTIAL_LIBEXECDIR}\")
    
    set(iSAAC_VERSION_FULL \"${iSAAC_VERSION_FULL}\")
    set(iSAAC_EXECUTABLE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
    set(iSAAC_LIBRARY_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
    ")

