#!/bin/bash
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
## file installBoost.sh
##
## Script to install boost
##
## author Come Raczy
##
################################################################################

REDIST_DIR=$1
INSTALL_DIR=$2
PARALLEL=$3
CXX=$4

. `dirname "$0"`/common.sh

BUILD_DIR=${INSTALL_DIR}/build
BIN_DIR=${INSTALL_DIR}/bin
LIB_DIR=${INSTALL_DIR}/lib
INCLUDE_DIR=${INSTALL_DIR}/include

SCRIPT=`basename "$0"`
VERSION=`echo ${iSAAC_BOOST_VERSION} | sed "s/\./_/g"`
SOURCE_TARBALL=${REDIST_DIR}/boost_${VERSION}.tar.bz2
TARBALL_COMPRESSION=j
SOURCE_DIR=${BUILD_DIR}/boost_${VERSION}

common_options $@

if [[ $CLEAN ]] ; then
    echo removing $SOURCE_DIR
    rm -rf $SOURCE_DIR ${INCLUDE_DIR}/boost ${LIB_DIR}/libboost_*.{a,so}
    exit 0
fi

common_create_source
#NO_BZIP2 disables dependency on libbz.h which is handy in cygwin
cd ${SOURCE_DIR} \
    && ./bootstrap.sh ${BOOTSTRAP_OPTIONS} --prefix=${INSTALL_DIR} --with-libraries=`echo ${iSAAC_BOOST_COMPONENTS} | sed "s/;/,/g"` \
    && echo "using gcc : : ${CXX} ;" >${SOURCE_DIR}/tools/build/v2/user-config.jam \
    && echo "modules.poke : NO_BZIP2 : 1 ;" >>${SOURCE_DIR}/tools/build/v2/user-config.jam \
    && ./bjam -q -j$PARALLEL ${BJAM_OPTIONS} --libdir=${INSTALL_DIR}/lib --layout=system link=static threading=multi install

if [ $? != 0 ] ; then echo "$SCRIPT: build failed: Terminating..." >&2 ; exit 1 ; fi

#echo "Cleaning up ${SOURCE_DIR}"  >&2
#rm -rf ${SOURCE_DIR}

echo "boost-$VERSION installed successfully"  >&2
