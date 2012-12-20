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
## file extractResults.sh
##
## extract meningful files from casava build folder into the location supplied
## as the second parameter
##
## author Roman Petrovski
##
################################################################################

#set -x
set -o pipefail
shopt -s compat31 2>/dev/null

LOGIFY_STDIN=@iSAAC_FULL_LIBEXECDIR@/logify.sh

CASAVA_BUILD_DIR=$1
shift

TARGET_RESULTS_DIR=$1
shift

#ISAAC_CP_VERBOSE=-v

TARGET_VARIANTS_DIR=${TARGET_RESULTS_DIR}/Parsed

for chrom in ${CASAVA_BUILD_DIR}/Parsed*/*; do

    chromName=$(basename $chrom)
    if [ "${chromName}" != "notMapped" ]; then
    
        chromTargetDir=${TARGET_VARIANTS_DIR}/$(basename $chrom)
        mkdir -p ${chromTargetDir} || exit 2
        find ${chrom} -maxdepth 1 -type f |xargs -I blah cp ${ISAAC_CP_VERBOSE} blah ${chromTargetDir}/ || exit 2
    
        chromIndelTargetDir=${chromTargetDir}/Indel
        [[ -d ${chromIndelTargetDir} ]] || mkdir ${chromIndelTargetDir} || exit 2
        [[ ! -d ${chrom}/Indel ]] || \
            find ${chrom}/Indel -maxdepth 1 -type f -not -name 'counts*' -not -name 'shadow*' -not -name '*fa.tmp' | \
                xargs -I blah cp ${ISAAC_CP_VERBOSE} blah ${chromIndelTargetDir}/ \
            || exit 2
        
        [[ ! -d ${chrom}/bam ]] || cp ${ISAAC_CP_VERBOSE} -r ${chrom}/bam ${chromTargetDir}/ || exit 2
        
        for sitesFile in ${chrom}/*/sites.txt.gz; do
            bin=$(basename $(dirname ${sitesFile}))
            mkdir -p ${chromTargetDir}/${bin} || exit 2
            cp ${ISAAC_CP_VERBOSE} ${sitesFile} ${chromTargetDir}/${bin}/ || exit 2
        done
        
    fi
done

cp ${ISAAC_CP_VERBOSE} -r ${CASAVA_BUILD_DIR}/html ${TARGET_RESULTS_DIR} || exit 2
cp ${ISAAC_CP_VERBOSE} -r ${CASAVA_BUILD_DIR}/stats ${TARGET_RESULTS_DIR} || exit 2
cp ${ISAAC_CP_VERBOSE} -r ${CASAVA_BUILD_DIR}/genome ${TARGET_RESULTS_DIR} || exit 2
