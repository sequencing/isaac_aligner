/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2012 Illumina, Inc.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Open Source Software License 1.
 **
 ** You should have received a copy of the Illumina Open Source
 ** Software License 1 along with this program. If not, see
 ** <https://github.com/downloads/sequencing/licenses/>.
 **
 ** The distribution includes the code libraries listed below in the
 ** 'redist' sub-directory. These are distributed according to the
 ** licensing terms governing each library.
 **
 ** \file isaac.cpp
 **
 ** \brief User-facing executable for alignment
 **
 ** \author Come Raczy
 **/
#include <sys/resource.h>

#include "common/Debug.hh"
#include "common/SystemCompatibility.hh"
#include "options/AlignOptions.hh"
#include "reference/ReferenceMetadata.hh"
#include "workflow/AlignWorkflowSerialization.hh"
#include "workflow/AlignWorkflow.hh"

void align(const isaac::options::AlignOptions &options);

int main(int argc, char *argv[])
{
    isaac::common::run(align, argc, argv);
}

void align(const isaac::options::AlignOptions &options)
{
    const unsigned long long availableMemory = options.memoryLimit * 1024 * 1024 * 1024;
    if (isaac::options::AlignOptions::memoryLimitUnlimited !=  options.memoryLimit)
    {
        ISAAC_THREAD_CERR << "align: Setting memory limit to " << availableMemory << " bytes." << std::endl;
        if (!isaac::common::ulimitV(availableMemory))
        {
            BOOST_THROW_EXCEPTION(isaac::common::ResourceException(
                errno, (boost::format("Failed to set the memory consumption limit to: %d bytes") % availableMemory).str()));
        }
    }

    isaac::workflow::AlignWorkflow workflow(
        options.argv,
        options.casavaArgv,
        options.flowcellLayoutList,
        options.barcodeMetadataList,
        options.allowVariableFastqLength,
        options.ignoreMissingBcls,
        options.ignoreMissingFilters,
        options.firstPassSeeds,
        0, //TODO: have a command-line argument to override the estimation-based value
        options.referenceMetadataList,
        options.tempDirectory,
        options.outputDirectory,
        options.jobs,
        options.repeatThreshold,
        options.mateDriftRange,
        options.neighborhoodSizeThreshold,
        availableMemory,
        options.ignoreNeighbors,
        options.ignoreRepeats,
        options.mapqThreshold,
        options.pfOnly,
        options.baseQualityCutoff,
        options.keepUnaligned,
        options.putUnalignedInTheBack,
        options.clipSemialigned,
        options.gappedMismatchesMax,
        options.scatterRepeats,
        options.gapMatchScore,
        options.gapMismatchScore,
        options.gapOpenScore,
        options.gapExtendScore,
        options.dodgyAlignmentScore,
        options.inputLoadersMax,
        options.tempSaversMax,
        options.tempLoadersMax,
        options.outputSaversMax,
        options.realignGaps,
        options.bamGzipLevel,
        options.expectedBgzfCompressionRatio,
        options.keepDuplicates,
        options.binRegexString,
        options.memoryControl,
        options.clusterIdList,
        options.userTemplateLengthStatistics,
        options.statsImageFormat);

    const boost::filesystem::path stateFilePath = options.tempDirectory / "AlignerState.xml";

    if (isaac::workflow::AlignWorkflow::Start != options.startFrom)
    {
        isaac::workflow::load(stateFilePath, workflow);
    }

    workflow.rewind(options.startFrom);
    isaac::workflow::AlignWorkflow::State targetState =
        (options.stopAt == isaac::workflow::AlignWorkflow::Last) ? workflow.getNextState() : options.stopAt;

    ISAAC_ASSERT_MSG(options.startFrom < targetState, "Target state must follow the start state");
    // save rewound state to avoid re-running on corrupt data if the next step breaks
    isaac::workflow::save(stateFilePath, workflow);

    while(targetState != workflow.step())
    {
        // save new state
        isaac::workflow::save(stateFilePath, workflow);
    }

    // save final state
    isaac::workflow::save(stateFilePath, workflow);
}

