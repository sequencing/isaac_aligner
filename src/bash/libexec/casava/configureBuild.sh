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
## <https://github.com/downloads/sequencing/licenses/>.
##
## The distribution includes the code libraries listed below in the
## 'redist' sub-directory. These are distributed according to the
## licensing terms governing each library.
##
################################################################################
##
## file configureBuild.sh
##
## Wrapper for launching CASAVA configureBuild.pl on iSAAC output data
##
## author Roman Petrovski
##
################################################################################

#set -x
set -o pipefail
shopt -s compat31 2>/dev/null

[[ "@iSAAC_CASAVA_BINDIR@/" == "/" ]] || iSAAC_CASAVA_BINDIR=@iSAAC_CASAVA_BINDIR@/

iSAAC_CASAVA_LIBEXECDIR=`${iSAAC_CASAVA_BINDIR}casavaConfig.sh --CASAVA_FULL_LIBEXECDIR` || exit $?

CASAVA_CONFIGURE_BUILD_PL_PATH=${iSAAC_CASAVA_BINDIR}configureBuild.pl
CASAVA_GENERATE_REPORTS_PL_PATH=${iSAAC_CASAVA_LIBEXECDIR}/generateReports.pl

argsCopy=("$@")

while (( ${#@} )); do
	param=$1
	shift
    if [[ $param == "--inSortedBam" || $param == "--ib" ]]; then
        inSortedBam=$1
        shift
    elif [[ $param == "--outDir" || $param == "--od" ]]; then
        outDir=$1
        shift
    elif [[ $param == "--help" || $param == "-h" ]]; then
        help=$param
    fi
done

# don't expect anything to happen if --help or -h is on the command line
[[ "" != "$help" ]] && $CASAVA_CONFIGURE_BUILD_PL_PATH $argsCopy && exit $?

[[ "" == "$outDir" || "" == "$inSortedBam" ]] && echo "ERROR: --inSortedBam and --outDir arguments are mandatory" >&2 && exit 2

 [[ ! -e "$inSortedBam" ]] && echo "ERROR: File not found: '$inSortedBam'" && exit 2
# [[ ! -e "${inSortedBam}.bai" ]] && \
#    echo -e "ERROR: File not found: '${inSortedBam}.bai'.\nPlease execute the following command:\nsamtools index $inSortedBam" && exit 2

bamDir=$(dirname "$inSortedBam")
matchSelectorStatsXml="$bamDir/../../../Stats/MatchSelectorStats.xml"
[[ ! -e "$matchSelectorStatsXml" ]] && echo "ERROR: File not found: '$matchSelectorStatsXml'" && exit 2

buildStatsXml="$bamDir/../../../Stats/BuildStats.xml"
[[ ! -e "$buildStatsXml" ]] && echo "ERROR: File not found: '$buildStatsXml'" && exit 2

ISAAC_SAMPLE_NAME=${ISAAC_SAMPLE_NAME-"$(basename $bamDir)"}
ISAAC_PROJECT_NAME=${ISAAC_PROJECT_NAME-"$(basename $(dirname $bamDir))"}

$CASAVA_CONFIGURE_BUILD_PL_PATH "${argsCopy[@]}" || exit $?

xsltproc \
	--stringparam ISAAC_SAMPLE_NAME_PARAM "$ISAAC_SAMPLE_NAME" \
	--stringparam ISAAC_PROJECT_NAME_PARAM "$ISAAC_PROJECT_NAME" \
	@iSAAC_FULL_DATADIR@/xsl/casava/MatchSelectorStatsToReadsIdx.xsl \
    $matchSelectorStatsXml >"$outDir/stats/Reads.idx" || exit $?

xsltproc \
	--stringparam ISAAC_SAMPLE_NAME_PARAM "$ISAAC_SAMPLE_NAME" \
	--stringparam ISAAC_PROJECT_NAME_PARAM "$ISAAC_PROJECT_NAME" \
	@iSAAC_FULL_DATADIR@/xsl/casava/BuildStatsToDupCountSummaryTxt.xsl \
    $buildStatsXml >"$outDir/stats/dupCount.summary.txt" || exit $?

# generateReports.pl tends to fail on bad data. Don't stop configuring because of that.
$CASAVA_GENERATE_REPORTS_PL_PATH --projectDir="$outDir" --target=duplicates || exit 0

