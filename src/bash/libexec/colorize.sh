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
