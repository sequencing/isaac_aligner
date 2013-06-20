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
## file check-source.sh
##
## Basic sanity checks on the source code before running the cppuint tests
##
## author Come Raczy
##
################################################################################

target=$1
shift
good=yes
for file in $* ; do
    check=`grep -nH CPPUNIT_TEST_SUITE_NAMED_REGISTRATION $file | grep -v registryName`
    if [[ $check ]] ; then
        if [[ $good ]] ; then
            good=
            echo >&2
            echo use of unchecked registry names: >&2
        fi
        echo "    "$check >&2
        echo >&2
    fi
done
if [[ $good ]] ; then
    echo checked > $target
else 
    exit 1
fi

