/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2012 Illumina, Inc.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Open Source Software License 1.
 **
 ** You should have received a copy of the Illumina Open Source
 ** Software License 1 along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** The distribution includes the code libraries listed below in the
 ** 'redist' sub-directory. These are distributed according to the
 ** licensing terms governing each library.
 **
 ** \file SelectMatchesTransition.cpp
 **
 ** Component to select the best matches among all possible candidates.
 **
 ** \author Come Raczy
 **/

#include <numeric>
#include <fstream>
#include <cerrno>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/thread.hpp>

#include "alignment/Alignment.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/FastIo.hh"
#include "common/ParallelSort.hpp"
#include "reference/Contig.hh"
#include "reference/ContigLoader.hh"

#include "alignment/matchSelector/MatchSelectorStatsXml.hh"

#include "workflow/alignWorkflow/SelectMatchesTransition.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

bool orderByTotalReadLengthDesc(const flowcell::FlowcellLayoutList &flowcellLayoutList,
                                const flowcell::TileMetadata &left, const flowcell::TileMetadata &right)
{
    const unsigned leftTotalReadLength = flowcell::getTotalReadLength(flowcellLayoutList.at(left.getFlowcellIndex()).getReadMetadataList());
    const unsigned rightTotalReadLength = flowcell::getTotalReadLength(flowcellLayoutList.at(right.getFlowcellIndex()).getReadMetadataList());
    // Also keep the natural order of tiles when the read lenghts are the same so that it is easier to track the
    // progress by monitoring the log output
    return  leftTotalReadLength > rightTotalReadLength ||
        (leftTotalReadLength == rightTotalReadLength && left.getIndex() < right.getIndex());
}

const flowcell::TileMetadataList sortByTotalReadLengthDesc(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    flowcell::TileMetadataList tileMetadataList)
{
    std::sort(tileMetadataList.begin(), tileMetadataList.end(),
              boost::bind(&orderByTotalReadLengthDesc, boost::ref(flowcellLayoutList), _1, _2));
    return tileMetadataList;
}


SelectMatchesTransition::SelectMatchesTransition(
        const unsigned ioOverlapParallelization,
        alignment::matchSelector::FragmentStorage &fragmentStorage,
        const alignment::MatchDistribution &matchDistribution,
        const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const boost::filesystem::path &tempDirectory,
        const unsigned int maxThreadCount,
        const TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned repeatThreshold,
        const int mateDriftRange,
        const bool allowVariableFastqLength,
        const bool cleanupIntermediary,
        const bool ignoreMissingBcls,
        const bool ignoreMissingFilters,
        const unsigned inputLoadersMax,
        const unsigned tempLoadersMax,
        const unsigned tempSaversMax,
        const alignment::MatchTally &matchTally,
        const alignment::TemplateLengthStatistics &userTemplateLengthStatistics,
        const unsigned mapqThreshold,
        const bool perTileTls,
        const bool pfOnly,
        const unsigned baseQualityCutoff,
        const bool keepUnaligned,
        const bool clipSemialigned,
        const bool clipOverlapping,
        const bool scatterRepeats,
        const unsigned gappedMismatchesMax,
        const bool avoidSmithWaterman,
        const int gapMatchScore,
        const int gapMismatchScore,
        const int gapOpenScore,
        const int gapExtendScore,
        const int minGapExtendScore,
        const unsigned semialignedGapLimit,
        const alignment::TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore,
        const bool qScoreBin,
        const boost::array<char, 256> &fullBclQScoreTable,
        const bool extractClusterXy
    )
    : matchLoadThreads_(tempLoadersMax),
      inputLoaderThreads_(inputLoadersMax),
      tileMetadataList_(tileMetadataList),
      processOrderTileMetadataList_(sortByTotalReadLengthDesc(flowcellLayoutList, tileMetadataList)),
      flowcellLayoutList_(flowcellLayoutList),
      ioOverlapThreads_(ioOverlapParallelization),
      nextUnprocessedTile_(processOrderTileMetadataList_.begin()),
      loadSlotAvailable_(true),
      flushSlotAvailable_(true),
      computeSlotAvailable_(true),

      matchTally_(matchTally),
      threadMatches_(ioOverlapParallelization, std::vector<alignment::Match>(getMaxTileMatches(matchTally_))),
      fragmentStorage_(fragmentStorage),
      matchLoader_(matchLoadThreads_),
      threadBclData_(ioOverlapParallelization, alignment::BclClusters(flowcell::getMaxTotalReadLength(flowcellLayoutList_) + flowcell::getMaxBarcodeLength(flowcellLayoutList_))),
      bclBaseCallsSource_(
          flowcellLayoutList_.end() == std::find_if(
              flowcellLayoutList_.begin(), flowcellLayoutList_.end(),
              boost::bind(&flowcell::Layout::getFormat, _1) == flowcell::Layout::Bcl) ? 0 :
                  new BclBaseCallsSource(
                      flowcellLayoutList_,
                      tileMetadataList_,
                      ignoreMissingBcls,
                      ignoreMissingFilters,
                      inputLoaderThreads_,
                      inputLoadersMax,
                      extractClusterXy)),
      fastqBaseCallsSource_(
          flowcellLayoutList_.end() == std::find_if(
              flowcellLayoutList_.begin(), flowcellLayoutList_.end(),
              boost::bind(&flowcell::Layout::getFormat, _1) == flowcell::Layout::Fastq) ? 0 :
                  new FastqBaseCallsSource(
                      flowcellLayoutList_,
                      tileMetadataList_,
                      allowVariableFastqLength,
                      inputLoaderThreads_,
                      inputLoadersMax)),
      bamBaseCallsSource_(
          flowcellLayoutList_.end() == std::find_if(
              flowcellLayoutList_.begin(), flowcellLayoutList_.end(),
              boost::bind(&flowcell::Layout::getFormat, _1) == flowcell::Layout::Bam) ? 0 :
                  new BamBaseCallsSource(
                      tempDirectory,
                      flowcellLayoutList_,
                      tileMetadataList_,
                      cleanupIntermediary,
                      inputLoaderThreads_,
                      inputLoadersMax)),
      bclBgzfBaseCallsSource_(
          flowcellLayoutList_.end() == std::find_if(
              flowcellLayoutList_.begin(), flowcellLayoutList_.end(),
              boost::bind(&flowcell::Layout::getFormat, _1) == flowcell::Layout::BclBgzf) ? 0 :
                  new BclBgzfBaseCallsSource(
                      flowcellLayoutList_,
                      tileMetadataList_,
                      ignoreMissingBcls,
                      ignoreMissingFilters,
                      inputLoaderThreads_,
                      inputLoadersMax,
                      extractClusterXy)),
      matchSelector_(
          fragmentStorage,
          matchDistribution,
          sortedReferenceMetadataList,
        maxThreadCount,
        tileMetadataList,
        barcodeMetadataList,
        flowcellLayoutList,
        repeatThreshold,
        mateDriftRange,
        userTemplateLengthStatistics,
        mapqThreshold,
        perTileTls,
        pfOnly,
        baseQualityCutoff,
        keepUnaligned,
        clipSemialigned,
        clipOverlapping,
        scatterRepeats,
        gappedMismatchesMax,
        avoidSmithWaterman,
        gapMatchScore,
        gapMismatchScore,
        gapOpenScore,
        gapExtendScore,
        minGapExtendScore,
        semialignedGapLimit,
        dodgyAlignmentScore),
        qScoreBin_(qScoreBin),
        fullBclQScoreTable_(fullBclQScoreTable),
        forceTermination_(false)
{
    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions constructor begin ")

    matchLoader_.reservePathBuffers(matchTally_.getMaxFilePathLength());

    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions before bclMapper_.reserveClusters ")
    BOOST_FOREACH(alignment::BclClusters &bclData, threadBclData_)
    {
        bclData.reserveClusters(flowcell::getMaxTileClusters(tileMetadataList_), extractClusterXy);
    }
    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions after bclData.reserveClusters ")

    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions constructor end ")
    ISAAC_THREAD_CERR << "Constructed the SelectMatchesTransition" << std::endl;
}

void SelectMatchesTransition::selectMatches(
    const common::ScoopedMallocBlock::Mode memoryControl,
    const boost::filesystem::path &matchSelectorStatsXmlPath,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics)
{
    {
        common::ScoopedMallocBlock  mallocBlock(memoryControl);
        nextUnprocessedTile_ = processOrderTileMetadataList_.begin();
        ioOverlapThreads_.execute(boost::bind(&SelectMatchesTransition::selectTileMatches, this, _1,
                                              boost::ref(matchTally_), boost::ref(barcodeTemplateLengthStatistics), boost::ref(mallocBlock)));

        ISAAC_ASSERT_MSG(loadSlotAvailable_ && computeSlotAvailable_ && flushSlotAvailable_,
                         "All slots must be available after the processing threads are gone");
        matchSelector_.unreserve();
    }

    matchSelector_.dumpStats(matchSelectorStatsXmlPath);


}


inline bool sortByTileBarcodeClusterLocation(const alignment::Match &lhs, const alignment::Match &rhs)
{
    return
        lhs.seedId.getTileBarcodeCluster() < rhs.seedId.getTileBarcodeCluster() ||
        (
            (lhs.seedId.getTileBarcodeCluster() == rhs.seedId.getTileBarcodeCluster()) && (lhs.location < rhs.location ||
                // same location for different seed numbers designates different locations for the fragment
                // ensure consistency of results by ordering by the seed as well
                (lhs.location == rhs.location && (lhs.seedId.getSeed() < rhs.seedId.getSeed())
                )
            )
        );
}

void SelectMatchesTransition::loadClusters(
    const unsigned threadNumber,
    const flowcell::TileMetadata &tileMetadata)
{
    const flowcell::Layout &flowcell = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex());
    if (flowcell::Layout::Fastq == flowcell.getFormat())
    {
        fastqBaseCallsSource_->loadClusters(tileMetadata, threadBclData_[threadNumber]);
    }
    else if (flowcell::Layout::Bam == flowcell.getFormat())
    {
        bamBaseCallsSource_->loadClusters(tileMetadata, threadBclData_[threadNumber]);
    }
    else if (flowcell::Layout::BclBgzf == flowcell.getFormat())
    {
        bclBgzfBaseCallsSource_->loadClusters(processOrderTileMetadataList_, tileMetadata, threadBclData_[threadNumber]);
    }
    else
    {
        ISAAC_ASSERT_MSG(flowcell::Layout::Bcl == flowcell.getFormat(), "Unsupported flowcell layout format " << flowcell.getFormat());
        bclBaseCallsSource_->loadClusters(tileMetadata, threadBclData_[threadNumber]);
    }

    // Bin QScore
    if(qScoreBin_)
    {
        ISAAC_THREAD_CERR << "Binning qscores" << std::endl;
        ISAAC_ASSERT_MSG(fullBclQScoreTable_.size() == 256, "QScore bin table incorrect size");
        alignment::BclClusters &bclData = threadBclData_[threadNumber];

        for(std::vector<char>::iterator itr = bclData.cluster(0); itr != bclData.end(); ++itr)
        {
            const int xx = int((unsigned char)*itr);
            *itr = fullBclQScoreTable_[xx];
        }
        ISAAC_THREAD_CERR << "Binning qscores done" << std::endl;
    }
}

void SelectMatchesTransition::selectTileMatches(
    const unsigned threadNumber,
    const alignment::MatchTally &matchTally,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    common::ScoopedMallocBlock &mallocBlock)
{
    while (true)
    {
        acquireLoadSlot();
        if (processOrderTileMetadataList_.end() == nextUnprocessedTile_)
        {
            releaseLoadSlot(false);
            return;
        }

        const flowcell::TileMetadata &tileMetadata = *nextUnprocessedTile_++;

        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&SelectMatchesTransition::releaseLoadSlot, this, _1))
        {

            ISAAC_THREAD_CERR << "Loading matches for " << tileMetadata << std::endl;
            matchLoader_.load(matchTally.getFileTallyList(tileMetadata), threadMatches_[threadNumber]);
            ISAAC_THREAD_CERR << "Loading matches done for " << tileMetadata << std::endl;

            if(threadMatches_[threadNumber].empty())
            {
                // The processing code below does not handle empty data too well.
                continue;
            }

            loadClusters(threadNumber, tileMetadata);
        }

            // sort the matches by SeedId and reference position.
//            ISAAC_THREAD_CERR << "Sorting matches for " << tileMetadata << std::endl;
//            std::sort(threadMatches_[threadNumber].begin(), threadMatches_[threadNumber].end(), sortByTileBarcodeClusterLocation);
//            ISAAC_THREAD_CERR << "Sorting matches done for " << tileMetadata << std::endl;

        acquireComputeSlot();
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&SelectMatchesTransition::releaseComputeSlot, this, _1))
        {
            // sort the matches by SeedId and reference position
            ISAAC_THREAD_CERR << "Sorting matches by barcode for " << tileMetadata << std::endl;
            {
                common::ScoopedMallocBlockUnblock unblock(mallocBlock);
                common::parallelSort(threadMatches_[threadNumber], sortByTileBarcodeClusterLocation);
            }
            ISAAC_THREAD_CERR << "Sorting matches by barcode done for " << tileMetadata << std::endl;

            matchSelector_.parallelSelect(matchTally, barcodeTemplateLengthStatistics, tileMetadata, threadMatches_[threadNumber], threadBclData_[threadNumber]);

            // There are only two sets of thread fragment dispatcher buffers (the one being flushed and the one we've just filled)
            // Wait for exclusive flush buffers access and swap the buffers before giving up the compute slot
            acquireFlushSlot();
            fragmentStorage_.prepareFlush();
        }

        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&SelectMatchesTransition::releaseFlushSlot, this, _1))
        {
            // now we can do out-of-sync flush while other thread does its compute
            fragmentStorage_.flush();
        }
    }
}


} // namespace alignWorkflow
} // namespace alignment
} // namespace isaac
