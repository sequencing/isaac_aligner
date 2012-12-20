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
        alignment::matchSelector::FragmentStorage &fragmentStorage,
        const reference::SortedReferenceXmlList &sortedReferenceXmlList,
        const bfs::path &tempDirectory,
        const unsigned int maxThreadCount,
        const TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned repeatThreshold,
        const int mateDriftRange,
        const bool allowVariableFastqLength,
        const bool ignoreMissingBcls,
        const bool ignoreMissingFilters,
        const unsigned inputLoadersMax,
        const unsigned tempLoadersMax,
        const unsigned tempSaversMax,
        const alignment::MatchTally &matchTally,
        const alignment::TemplateLengthStatistics &userTemplateLengthStatistics,
        const unsigned mapqThreshold,
        const bool pfOnly,
        const unsigned baseQualityCutoff,
        const bool keepUnaligned,
        const bool clipSemialigned,
        const unsigned gappedMismatchesMax,
        const bool scatterRepeats,
        const int gapMatchScore,
        const int gapMismatchScore,
        const int gapOpenScore,
        const int gapExtendScore,
        const alignment::TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore
    )
    : bclLoadThreads_(inputLoadersMax),
      matchLoadThreads_(tempLoadersMax),
      tileMetadataList_(tileMetadataList),
      processOrderTileMetadataList_(sortByTotalReadLengthDesc(flowcellLayoutList, tileMetadataList)),
      flowcellLayoutList_(flowcellLayoutList),
      ioOverlapThreads_(processingStages_),
      nextUnprocessedTile_(processOrderTileMetadataList_.begin()),
      loadSlotAvailable_(true),
      flushSlotAvailable_(true),
      computeSlotAvailable_(true),

      matchTally_(matchTally),
      threadMatches_(processingStages_, std::vector<alignment::Match>(getMaxTileMatches(matchTally_))),
      fragmentStorage_(fragmentStorage),
      bclMapper_(ignoreMissingBcls, flowcell::getMaxTotalReadLength(flowcellLayoutList_), bclLoadThreads_, flowcell::getMaxTileClulsters(tileMetadataList_)),
      fastqLoader_(allowVariableFastqLength),
      matchLoader_(matchLoadThreads_),
      threadFastqFilePaths_(processingStages_, std::vector<boost::filesystem::path>(2)), //read 1 and read 2 paths
      threadBclFilePaths_(processingStages_, std::vector<boost::filesystem::path>(flowcell::getMaxTotalReadLength(flowcellLayoutList_))),
      threadBclData_(processingStages_, alignment::BclClusters(flowcell::getMaxTotalReadLength(flowcellLayoutList_))),
      threadFilterFilePaths_(processingStages_),
      filtersMapper_(ignoreMissingFilters),
      matchSelector_(fragmentStorage,
        sortedReferenceXmlList,
        maxThreadCount,
        tileMetadataList,
        barcodeMetadataList,
        flowcellLayoutList,
        repeatThreshold,
        mateDriftRange,
        userTemplateLengthStatistics,
        mapqThreshold,
        pfOnly,
        baseQualityCutoff,
        keepUnaligned,
        clipSemialigned,
        gappedMismatchesMax,
        scatterRepeats,
        gapMatchScore,
        gapMismatchScore,
        gapOpenScore,
        gapExtendScore,
        dodgyAlignmentScore)
{
    const boost::filesystem::path longestBaseCallsPath =  flowcell::getLongestBaseCallsPath(tileMetadataList_);
    const unsigned highestTileNumber = std::max_element(tileMetadataList_.begin(), tileMetadataList_.end(),
                                                        boost::bind(&flowcell::TileMetadata::getTile, _1)<
                                                            boost::bind(&flowcell::TileMetadata::getTile, _2))->getTile();

    const bool compressedFound =
        tileMetadataList_.end() != std::find_if(tileMetadataList_.begin(), tileMetadataList_.end(),
                                                boost::bind(&flowcell::TileMetadata::getCompression, _1) !=
                                                    flowcell::TileMetadata::NoCompression);

    // reserve memory needed for fastq processing
    boost::filesystem::path longestFastqFilePath;
    flowcell::Layout::getFastqFilePath(1, 1, longestBaseCallsPath, flowcell::TileMetadata::GzCompression, longestFastqFilePath);
    BOOST_FOREACH(std::vector<boost::filesystem::path> &vp, threadFastqFilePaths_)
    {
        BOOST_FOREACH(boost::filesystem::path &p, vp)
        {
            // this has to be done separately for each path or else they all share one buffer
            p = longestFastqFilePath.c_str();
        }
    }

    fastqLoader_.reservePathBuffers(longestFastqFilePath.string().size());
    // reserve memory needed for bcl processing
    const unsigned insanelyHighCycleNumber = 9999;
    boost::filesystem::path longestBclFilePath;
    flowcell::Layout::getBclFilePath(highestTileNumber, 1, longestBaseCallsPath, insanelyHighCycleNumber,
                                     flowcell::TileMetadata::GzCompression, longestBclFilePath);

    bclMapper_.reserveBuffers(longestBclFilePath.string().size(), compressedFound);


    matchLoader_.reservePathBuffers(alignment::MatchTally::getMaxFilePathLength(tempDirectory));
    BOOST_FOREACH(std::vector<boost::filesystem::path> &vp, threadBclFilePaths_)
    {
        BOOST_FOREACH(boost::filesystem::path &p, vp)
        {
            // this has to be done separately for each path or else they all share one buffer
            p = longestBclFilePath.c_str();
        }
    }

    BOOST_FOREACH(alignment::BclClusters &bclData, threadBclData_)
    {
        bclData.reserveClusters(flowcell::getMaxTileClulsters(tileMetadataList_));
    }

    boost::filesystem::path longestFilterFilePath = flowcell::Layout::getLongestFilterFilePath(flowcellLayoutList_);
    BOOST_FOREACH(boost::filesystem::path &p, threadFilterFilePaths_)
    {
        // this has to be done separately for each path or else they all share one buffer
        p = longestFilterFilePath.c_str();
    }
    filtersMapper_.reservePathBuffers(longestFilterFilePath.string().size());
    filtersMapper_.reserveBuffer(flowcell::getMaxTileClulsters(tileMetadataList_));

    ISAAC_THREAD_CERR << "Constructed the match selector" << std::endl;
}

void SelectMatchesTransition::selectMatches(
    const common::ScoopedMallocBlock::Mode memoryControl,
    const bfs::path &matchSelectorStatsXmlPath)
{
    {
        common::ScoopedMallocBlock  mallocBlock(memoryControl);
        nextUnprocessedTile_ = processOrderTileMetadataList_.begin();
        ioOverlapThreads_.execute(boost::bind(&SelectMatchesTransition::selectTileMatches, this, _1,
                                              boost::ref(matchTally_), boost::ref(mallocBlock)));

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
    if (flowcell::Layout::Fastq == flowcell.getFormat() || flowcell::Layout::FastqGz == flowcell.getFormat())
    {
        ISAAC_THREAD_CERR << "Loading Fastq data for " << tileMetadata << std::endl;

        flowcell.getFastqFilePath(1, tileMetadata.getLane(), tileMetadata.getBaseCallsPath(),
                                  flowcell::Layout::FastqGz == flowcell.getFormat(),
                                  threadFastqFilePaths_[threadNumber].at(0));
        flowcell.getFastqFilePath(2, tileMetadata.getLane(), tileMetadata.getBaseCallsPath(),
                                  flowcell::Layout::FastqGz == flowcell.getFormat(),
                                  threadFastqFilePaths_[threadNumber].at(1));
        fastqLoader_.open(threadFastqFilePaths_[threadNumber][0], threadFastqFilePaths_[threadNumber][1]);

        alignment::BclClusters &bclData = threadBclData_[threadNumber];

        const unsigned clustersToLoad = tileMetadata.getClusterCount();
        ISAAC_THREAD_CERR << "Resetting Fastq data for " << clustersToLoad << " clusters" << std::endl;
        bclData.reset(flowcell::getTotalReadLength(flowcell.getReadMetadataList()), clustersToLoad);
        ISAAC_THREAD_CERR << "Resetting Fastq data done for " << bclData.getClusterCount() << " clusters" << std::endl;

        std::vector<char>::iterator clustersEnd = bclData.cluster(0);
        unsigned clustersLoaded = fastqLoader_.loadClusters(
            bclData.getClusterCount(), flowcell.getReadMetadataList(), clustersEnd);
        ISAAC_ASSERT_MSG(clustersToLoad == clustersLoaded, "Loaded mismatching number of clusters: " << clustersLoaded << " expected: " << tileMetadata);

        bclData.pf().clear();
        for(unsigned cluster = 0; clustersLoaded > cluster; ++cluster)
        {
            bclData.pf().push_back(true);
        }


        ISAAC_THREAD_CERR << "Loading Fastq data done. Loaded " << clustersLoaded << " clusters for " << tileMetadata << std::endl;

    }
    else
    {
        ISAAC_THREAD_CERR << "Loading Bcl data for " << tileMetadata << std::endl;

        const std::vector<unsigned> &cycleList = flowcell.getAllCycleNumbers();
        ISAAC_ASSERT_MSG(threadBclFilePaths_[threadNumber].size() >= cycleList.size(), "tiles expected to be ordered so that number of cycles goes down")
        threadBclFilePaths_[threadNumber].resize(cycleList.size());
        unsigned pos = 0;
        BOOST_FOREACH(const unsigned int cycle, cycleList)
        {
            flowcell::Layout::getBclFilePath(
                tileMetadata.getTile(), tileMetadata.getLane(),
                tileMetadata.getBaseCallsPath(), cycle,
                tileMetadata.getCompression(), threadBclFilePaths_[threadNumber][pos++]);
        }

        bclMapper_.mapTile(threadBclFilePaths_[threadNumber], tileMetadata.getClusterCount());
        ISAAC_THREAD_CERR << "Loading Bcl data done for " << tileMetadata << std::endl;

        ISAAC_THREAD_CERR << "Loading Filter data for " << tileMetadata << std::endl;
        flowcell.getFiltersFilePath(tileMetadata.getTile(), tileMetadata.getLane(), tileMetadata.getBaseCallsPath(),
                                    threadFilterFilePaths_[threadNumber]);
        filtersMapper_.mapTile(threadFilterFilePaths_[threadNumber], tileMetadata.getClusterCount());
        ISAAC_THREAD_CERR << "Loading Filter data done for " << tileMetadata << std::endl;

        // bclToClusters mainly does transposition of bcl cycles to clusters which is a non-io operation.
        // However, the amount of CPU required is relatively low, and occurs on a single thread.
        // Avoid locking all the cores for the duration of this...
        // Also, bclMapper_ and filtersMapper_ are shared between the threads at the moment.
        bclToClusters(bclMapper_, filtersMapper_, tileMetadata, threadBclData_[threadNumber]);
    }
}

void SelectMatchesTransition::selectTileMatches(const unsigned threadNumber,
                                      const alignment::MatchTally &matchTally,
                                      common::ScoopedMallocBlock &mallocBlock)
{
    while (true)
    {
        // Note!!!, the concurrency locking is not exception-safe. The processing is exepected to terminate
        // in case of exception escaping from this method!
        acquireLoadSlot();

        if (processOrderTileMetadataList_.end() == nextUnprocessedTile_)
        {
            releaseLoadSlot();
            return;
        }

        const flowcell::TileMetadata &tileMetadata = *nextUnprocessedTile_++;

        // this scope controls the memory consumption for the loaded tile bcl and match data
        {
            ISAAC_THREAD_CERR << "Loading matches for " << tileMetadata << std::endl;
            matchLoader_.load(matchTally.getFileTallyList(tileMetadata), threadMatches_[threadNumber]);
            ISAAC_THREAD_CERR << "Loading matches done for " << tileMetadata << std::endl;

            if(threadMatches_[threadNumber].empty())
            {
                // The processing code below does not handle empty data too well.
                releaseLoadSlot();
                continue;
            }

            loadClusters(threadNumber, tileMetadata);

            releaseLoadSlot();

            // sort the matches by SeedId and reference position.
//            ISAAC_THREAD_CERR << "Sorting matches for " << tileMetadata << std::endl;
//            std::sort(threadMatches_[threadNumber].begin(), threadMatches_[threadNumber].end(), sortByTileBarcodeClusterLocation);
//            ISAAC_THREAD_CERR << "Sorting matches done for " << tileMetadata << std::endl;

            acquireComputeSlot();
            // sort the matches by SeedId and reference position
            ISAAC_THREAD_CERR << "Sorting matches by barcode for " << tileMetadata << std::endl;
            {
                common::ScoopedMallocBlockUnblock unblock(mallocBlock);
                common::parallelSort(threadMatches_[threadNumber], sortByTileBarcodeClusterLocation);
            }
            ISAAC_THREAD_CERR << "Sorting matches by barcode done for " << tileMetadata << std::endl;

            matchSelector_.parallelSelect(matchTally, tileMetadata, threadMatches_[threadNumber], threadBclData_[threadNumber]);
        }

        // There are only two sets of thread fragment dispatcher buffers (the one being flushed and the one we've just filled)
        // Wait for exclusive flush buffers access and swap the buffers before giving up the compute slot
        acquireFlushSlot();
        fragmentStorage_.prepareFlush();

        releaseComputeSlot();

        // now we can do out-of-sync flush while other thread does its compute
        fragmentStorage_.flush();
        releaseFlushSlot();
    }
}

void SelectMatchesTransition::bclToClusters(
    const io::ParallelBclMapper &bclMapper,
    const io::FiltersMapper &filtersMapper,
    const flowcell::TileMetadata &tileMetadata,
    alignment::BclClusters &bclData) const
{

    ISAAC_THREAD_CERR << "Resetting Bcl data for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
    bclData.reset(bclMapper.getCyclesCount(), tileMetadata.getClusterCount());
    ISAAC_THREAD_CERR << "Resetting Bcl data done for " << bclData.getClusterCount() << " bcl clusters" << std::endl;


    ISAAC_THREAD_CERR << "Transposing Bcl data for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
    const clock_t startTranspose = clock();
    bclMapper.transpose(bclData.cluster(0));
    ISAAC_THREAD_CERR << "Transposing Bcl data done for " << bclData.getClusterCount() << " bcl clusters in " << (clock() - startTranspose) / 1000 << "ms" << std::endl;

    ISAAC_THREAD_CERR << "Extracting Pf values for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
    bclData.pf().clear();
    filtersMapper.getPf(std::back_inserter(bclData.pf()));
    assert(bclData.pf().size() == tileMetadata.getClusterCount());
    ISAAC_THREAD_CERR << "Extracting Pf values done for " << bclData.getClusterCount() << " bcl clusters" << std::endl;
}

} // namespace alignWorkflow
} // namespace alignment
} // namespace isaac
