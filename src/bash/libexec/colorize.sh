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
## file colorize.sh
##
## stream colorizer
##
## author Roman Petrovski
##
################################################################################

#set -x
set -o pipefail
shopt -s compat31 2>/dev/null

EXPRESSION=$1

while read l; do 
	(echo $l | grep -i -E "$EXPRESSION" --color || echo $l)
done
