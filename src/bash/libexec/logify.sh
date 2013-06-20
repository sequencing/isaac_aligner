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
## file logifystdin.sh
##
## logging stream formatter
##
## author Roman Petrovski
##
################################################################################

#set -x
set -o pipefail
shopt -s compat31 2>/dev/null

LOG_LEVEL=$1
shift

LOG_DATE_FORMAT=$1
shift

TARGET=$1
shift

case $LOG_LEVEL in
    '' )
        cat ;;
	0 )
		cat >/dev/null;;
	1 )
		(grep --line-buffered -i -E '(info|warning|error)' || [[ 1 == $? ]]) | while read -r l; do
            dt=$(date "+$LOG_DATE_FORMAT") || exit 2
			(echo -en [$dt]"\t[$HOSTNAME]\t[$TARGET]\t" && echo "$l") || exit 2
		done;;
	* )
		while read -r l; do
			dt=$(date "+$LOG_DATE_FORMAT") || exit 2
			(echo -en [$dt]"\t[$HOSTNAME]\t[$TARGET]\t" && echo "$l") || exit 2
		done;;
esac

