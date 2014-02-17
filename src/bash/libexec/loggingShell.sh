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
## file loggingshell.sh
##
## make $(SHELL) wrapper that passess the STDERR of a command line through the 
## logging stream formatter
##
## author Roman Petrovski
##
################################################################################

#set -x
set -o pipefail
shopt -s compat31 2>/dev/null

LOGIFY_STDIN=@iSAAC_HOME@@iSAAC_FULL_LIBEXECDIR@/logify.sh

# NFS default caching is 60 seconds + "small" amount of seconds for isilon
# The downside of long waits is that if there are some incorrectly-named phony
# targets in the workflow, we will wait long for each one.
# Not that we're supposed to have the incorrectly-named phony targets in our workflows...
iSAAC_MISSING_FILE_WAIT_DELAY=${iSAAC_MISSING_FILE_WAIT_DELAY:-75}

LOG_LEVEL=$1
shift

LOG_DATE_FORMAT=$1
shift

MAKE=$1
shift

FULL_TARGET=$1
#SHORT_TARGET=$(basename "$FULL_TARGET")
#[[ '.sentinel' == "$SHORT_TARGET" ]] && SHORT_TARGET=$(basename $(basename "$FULL_TARGET"))/.sentinel

# padding is needed so that if the FULL_TARGET is shorter than 30, it still gets properly extracted
padded_full_target="                                                            $FULL_TARGET"
SHORT_TARGET=${padded_full_target:(-50)}
SHORT_TARGET=${SHORT_TARGET// }
shift

ALL_PREREQS=$1
shift

# some distributed file sytems don't guarantee the immediate visibility of a file on a different node than the one
# that produces it. Try to sleep and generate warning if that's still the case. If your phony target generates this
# warning, append .done to the end of the target name.
for prereq in $ALL_PREREQS; do
    [[ ! "$prereq" =~ '\.done$|^all$|^clean$|^default$' ]] && [[ ! -e "$prereq" ]] && sleep ${iSAAC_MISSING_FILE_WAIT_DELAY} && [[ ! -e "$prereq" ]] && \
        echo WARNING: prereq "'$prereq'" is not in the file system after waiting for ${iSAAC_MISSING_FILE_WAIT_DELAY} seconds \
        | ${LOGIFY_STDIN} "$LOG_LEVEL" "${LOG_DATE_FORMAT}" "$SHORT_TARGET" >&2
done

NEWER_PREREQS=$1
shift

SHELL_ORI=$1
shift

[[ "" == "$FULL_TARGET" ]] && echo "ERROR: target cannot be emtpy. Reason: $NEWER_PREREQS\nDependencies: $ALL_PREREQS\nCmd: $@" >&2 && exit 2 

# get rid of -c that comes from make
[[ '-c' == $1 ]] && shift

removed_multilines=$(echo "$@" | while read -r l; do line_without_newline=$(echo -n "$l"); echo -n "${line_without_newline%\\} "; done)

(echo -en "$MAKE Target:\t$FULL_TARGET\n$MAKE Reason:\t$NEWER_PREREQS\n$MAKE Prereqs:\t$ALL_PREREQS\n$MAKE Cmd:\t" && echo "$removed_multilines") \
    | ${LOGIFY_STDIN} "$LOG_LEVEL" "${LOG_DATE_FORMAT}" "$SHORT_TARGET" >&2

#swap own stderr and stdout as child shell will have them swapped as well
exec 3>&1 1>&2 2>&3 3>&-

export iSAAC_LOG_LEVEL=$LOG_LEVEL;

3>&1 1>&2 2>&3 3>&- $SHELL_ORI -c "$@" | ${LOGIFY_STDIN} "$LOG_LEVEL" "${LOG_DATE_FORMAT}" "$SHORT_TARGET" 

exit $?
