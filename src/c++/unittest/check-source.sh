#!/bin/bash
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

