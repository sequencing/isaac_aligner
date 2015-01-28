/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** BSD 2-Clause License
 **
 ** You should have received a copy of the BSD 2-Clause License
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** \file FindMatchesTransition.cpp
 **
 ** \brief see FindMatchesTransition.hh
 **
 ** \author Roman Petrovski
 **/

#include "alignment/MatchFinder.hh"
#include "alignment/SeedLoader.hh"
#include "alignment/SeedMemoryManager.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/ParallelSort.hpp"
#include "demultiplexing/DemultiplexingStatsXml.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "workflow/alignWorkflow/BamDataSource.hh"
#include "workflow/alignWorkflow/BclBgzfDataSource.hh"
#include "workflow/alignWorkflow/BclDataSource.hh"
#include "workflow/alignWorkflow/FastqDataSource.hh"
#include "workflow/alignWorkflow/FindMatchesTransition.hh"


namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

FindMatchesTransition::FindMatchesTransition(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const bool allowVariableFastqLength,
    const bool cleanupIntermediary,
    const bool ignoreMissingBcls,
    const unsigned firstPassSeeds,
    const unsigned long availableMemory,
    const unsigned clustersAtATimeMax,
    const bfs::path &tempDirectory,
    const bfs::path &demultiplexingStatsXmlPath,
    const unsigned int maxThreadCount,
    const unsigned repeatThreshold,
    const unsigned neighborhoodSizeThreshold,
    const bool ignoreNeighbors,
    const bool ignoreRepeats,
    const unsigned inputLoadersMax,
    const unsigned tempSaversMax,
    const common::ScoopedMallocBlock::Mode memoryControl,
    const std::vector<size_t> &clusterIdList,
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList
    )
    : flowcellLayoutList_(flowcellLayoutList)
    , tempDirectory_(tempDirectory)
    , demultiplexingStatsXmlPath_(demultiplexingStatsXmlPath)
    , coresMax_(maxThreadCount)
    , repeatThreshold_(repeatThreshold)
    , neighborhoodSizeThreshold_(neighborhoodSizeThreshold)
    , barcodeMetadataList_(barcodeMetadataList)
    , allowVariableFastqLength_(allowVariableFastqLength)
    , cleanupIntermediary_(cleanupIntermediary)
    , ignoreMissingBcls_(ignoreMissingBcls)
    , firstPassSeeds_(firstPassSeeds)
    , availableMemory_(availableMemory)
    , clustersAtATimeMax_(clustersAtATimeMax)
    , ignoreNeighbors_(ignoreNeighbors)
    , ignoreRepeats_(ignoreRepeats)
    , inputLoadersMax_(inputLoadersMax)
    , tempSaversMax_(tempSaversMax)
    , memoryControl_(memoryControl)
    , clusterIdList_(clusterIdList)
    , sortedReferenceMetadataList_(sortedReferenceMetadataList)
    // Have thread pool for the maximum number of threads we may potentially need.
    , threads_(std::max(inputLoadersMax_, std::max(coresMax_, tempSaversMax_)))
{
}

std::vector<std::vector<unsigned> > FindMatchesTransition::getSeedIndexListPerIteration(const flowcell::Layout &flowcell) const
{
    std::vector<std::vector<unsigned> >seedIndexListPerIteration(maxIterations_); // assume two iterations
    std::vector<unsigned> countsPerRead(flowcell.getReadMetadataList().size(), 0);
    BOOST_FOREACH(const alignment::SeedMetadata &seedMetadata, flowcell.getSeedMetadataList())
    {
        assert(countsPerRead.size() > seedMetadata.getReadIndex());
        // allocate the first seed of each read to the first iteration, and all
        // subsequent seeds of the read to the second iteration
        const unsigned iteration = (firstPassSeeds_ > countsPerRead[seedMetadata.getReadIndex()]) ? 0 : 1;
        seedIndexListPerIteration[iteration].push_back(seedMetadata.getIndex());
        ++countsPerRead[seedMetadata.getReadIndex()];
    }
    // trim iterations without any seeds
    while (!seedIndexListPerIteration.empty() && seedIndexListPerIteration.back().empty())
    {
        seedIndexListPerIteration.resize(seedIndexListPerIteration.size() - 1);
    }
    return seedIndexListPerIteration;
}
static const unsigned standardOpenFileHandlesCount(3); // cin, cout, cerr

void FindMatchesTransition::resolveBarcodes(
    const flowcell::Layout &flowcell,
    const flowcell::BarcodeMetadataList &barcodeGroup,
    // this needs to contain only the tiles that we are processing by placed at the tile.getIndex()
    const flowcell::TileMetadataList &allTiles,
    BarcodeSource &barcodeSource,
    // this contains tiles we are processing but they are not placed at the tile.getIndex()
    flowcell::TileMetadataList unprocessedTiles,
    alignment::matchFinder::TileClusterInfo &tileClusterInfo,
    demultiplexing::DemultiplexingStats &demultiplexingStats)
{
    ISAAC_ASSERT_MSG(!barcodeGroup.empty(), "At least 'none' barcode must be defined");
    if (1 == barcodeGroup.size())
    {
        const flowcell::BarcodeMetadata barcode = barcodeGroup.at(0);
        ISAAC_ASSERT_MSG(barcode.isNoIndex(), "If barcode group has only one entry it must be the 'NoIndex' barcode");
        BOOST_FOREACH(const flowcell::TileMetadata &tile, unprocessedTiles)
        {
            for (unsigned clusterId = 0; clusterId < tile.getClusterCount(); ++clusterId)
            {
                tileClusterInfo.setBarcodeIndex(tile.getIndex(), clusterId, barcode.getIndex());
                demultiplexingStats.recordBarcode(demultiplexing::BarcodeId(tile.getIndex(), barcode.getIndex(), clusterId, 0));
            }
            ISAAC_THREAD_CERR << "Forced barcode index for clusters of " << tile << " to " << barcode << std::endl;
        }
    }
    else
    {
        demultiplexing::BarcodeResolver barcodeResolver(allTiles, barcodeMetadataList_, barcodeGroup);

        flowcell::TileMetadataList currentTiles; currentTiles.reserve(unprocessedTiles.size());

        while (!unprocessedTiles.empty())
        {
            currentTiles.clear();
            if (!demultiplexing::BarcodeMemoryManager::selectTiles(unprocessedTiles, currentTiles))
            {
                BOOST_THROW_EXCEPTION(common::MemoryException("Insufficient memory to load barcodes even for just one tile: " +
                    boost::lexical_cast<std::string>(unprocessedTiles.back())));
            }

            std::vector<demultiplexing::Barcode> barcodes;
            // this will take at most the same amount of ram as a set of singleseeds
            ISAAC_ASSERT_MSG(barcodeGroup.size(), "Barcode list must be not empty");
            ISAAC_ASSERT_MSG(barcodeGroup.at(0).isDefault(), "The very first barcode must be the 'unknown indexes or no index' one");
            barcodeSource.loadBarcodes(barcodeGroup.at(0).getIndex(), currentTiles, barcodes);
            barcodeResolver.resolve(barcodes, demultiplexingStats);

            BOOST_FOREACH(const demultiplexing::Barcode &barcode, barcodes)
            {
                tileClusterInfo.setBarcodeIndex(barcode.getTile(), barcode.getCluster(), barcode.getBarcode());
            }
        }
    }
}

/**
 * \brief Executes matchFinder in exact match mode.
 *
 * \param unprocessedTiles [inout] Processed tiles are removed from the list upon return
 * \param tileClusterInfo  [out]   Match finder marks read as complete if exact match is found.
 *                                 Barcode resolution set the barcode index for each cluster of the returned
 *                                 tile set.
 * \param foundMatches     [out]   Updated upon return.
 *
 * \return Returns the list of tiles for which the match finding was performed
 */
template <typename KmerT>
flowcell::TileMetadataList FindMatchesTransition::findSingleSeedMatches(
    const flowcell::Layout &flowcell,
    const std::vector<unsigned> &seedIndexList,
    const bool finalPass,
    flowcell::TileMetadataList &unprocessedTiles,
    alignment::matchFinder::TileClusterInfo &tileClusterInfo,
    SeedSource<KmerT> &seedSource,
    demultiplexing::DemultiplexingStats &demultiplexingStats,
    FoundMatchesMetadata &foundMatches)
{
    const unsigned seedLoaderOpenFileHandlesCount(inputLoadersMax_);

    alignment::SeedMetadataList seedMetadataList;
    BOOST_FOREACH(size_t index, seedIndexList)
    {
        seedMetadataList.push_back(flowcell.getSeedMetadataList().at(index));
    }
    ISAAC_THREAD_CERR << "Finding Single-seed matches for " << seedMetadataList << "with repeat threshold: " <<
        repeatThreshold_ << std::endl;

    alignment::MatchFinder<KmerT> matchFinder(sortedReferenceMetadataList_, tempDirectory_,
                            unprocessedTiles, flowcell.getReadMetadataList(), flowcell.getSeedMetadataList(),
                            0,
                            ignoreNeighbors_, ignoreRepeats_,
                            repeatThreshold_, neighborhoodSizeThreshold_,
                            foundMatches.matchTally_, tileClusterInfo, threads_, coresMax_, tempSaversMax_,
                            standardOpenFileHandlesCount + seedLoaderOpenFileHandlesCount);

    flowcell::TileMetadataList currentTiles; currentTiles.reserve(unprocessedTiles.size());

    seedSource.initBuffers(unprocessedTiles, seedMetadataList);

    alignment::SeedMemoryManager<KmerT> seedMemoryManager(
        barcodeMetadataList_, flowcell.getReadMetadataList(), seedMetadataList, unprocessedTiles);

    if (!seedMemoryManager.selectTiles(
        unprocessedTiles, tileClusterInfo, matchFinder.getMaxTileCount(), tempSaversMax_, currentTiles))
    {
        BOOST_THROW_EXCEPTION(common::MemoryException("Insufficient memory to load seeds even for just one tile: " +
            boost::lexical_cast<std::string>(unprocessedTiles.back())));
    }

    {
        std::vector<alignment::Seed<KmerT> > seeds;
        seedMemoryManager.allocate(currentTiles, seeds);

        common::ScoopedMallocBlock  mallocBlock(memoryControl_);
        seedSource.generateSeeds(currentTiles, tileClusterInfo, seeds, mallocBlock);
        matchFinder.setTiles(currentTiles);

        ISAAC_THREAD_CERR << "Finding Exact single-seed matches for " << seedMetadataList << "with repeat threshold: " <<
            repeatThreshold_ << std::endl;
        foundMatches.matchDistribution_.consolidate(
            matchFinder.findMatches(seeds.begin(), seedSource.getReferenceSeedBounds(), false, finalPass));
        ISAAC_THREAD_CERR << "Finding Exact single-seed matches done for " << seedMetadataList << std::endl;

        std::vector<alignment::Seed<KmerT> >().swap(seeds);
    }

    ISAAC_THREAD_CERR << "Finding Single-seed matches done for " << seedMetadataList << std::endl;

    return currentTiles;
}

/**
 * \brief Turns seeds belonging to complete cluster reads into N-seeds, sorts seeds so that
 *        N-seeds get moved to the back of each range.
 *
 * \param clusterInfo         information about which cluster reads can be masked
 */
template <typename KmerT>
void FindMatchesTransition::maskCompleteReadSeeds(
    const alignment::SeedMetadataList &allSeedMetadataList,
    const alignment::matchFinder::TileClusterInfo &clusterInfo,
    std::vector<alignment::Seed<KmerT> > &seeds,
    const std::vector<typename std::vector<alignment::Seed<KmerT> >::iterator> &referenceSeedBounds,
    common::ScoopedMallocBlock  &mallocBlock)
{
    typedef alignment::Seed<KmerT> SeedT;
    typedef typename std::vector<SeedT>::iterator SeedIterator;
    // the Seeds buffer might turn out to be bigger than what we need if some clusters map to barcodes which
    // have unmapped references. This is not perceived to be the major scenario, so, some unused memory is acceptable
    ISAAC_ASSERT_MSG(seeds.end() >= referenceSeedBounds.back(), "Computed end is past the end of the reserved buffer");

    SeedIterator referenceSeedsBegin = seeds.begin();
    BOOST_FOREACH(SeedIterator referenceSeedsEnd, referenceSeedBounds)
    {
        ISAAC_THREAD_CERR << "Masking " << std::distance(referenceSeedsBegin, referenceSeedsEnd) << " seeds" << std::endl;
        const clock_t startMask = clock();
        unsigned long masked = 0;
        BOOST_FOREACH(SeedT &seed, std::make_pair(referenceSeedsBegin, referenceSeedsEnd))
        {
            if (!seed.isNSeed() &&
                clusterInfo.isReadComplete(seed.getTile(), seed.getCluster(),
                                           allSeedMetadataList.at(seed.getSeedIndex()).getReadIndex()))
            {
                // dont' set them to be lowest n-seeds as lowest n-seeds get no-matches stored for them
                seed.makeNSeed(false);
                ++masked;
            }
        }
        ISAAC_THREAD_CERR << "Masking " << std::distance(referenceSeedsBegin, referenceSeedsEnd) <<
            " seeds done (" << masked <<
            " masked) in " << (clock() - startMask) / 1000 << "ms" << std::endl;

        ISAAC_THREAD_CERR << "Sorting " << std::distance(referenceSeedsBegin, referenceSeedsEnd) << " seeds" << std::endl;
        const clock_t startSort = clock();
        {
            common::ScoopedMallocBlockUnblock unblock(mallocBlock);
            // comparing the full kmer is required to push the N-seeds off to the very end.
            common::parallelSort(referenceSeedsBegin, referenceSeedsEnd, &alignment::orderByKmerSeedIndex<KmerT>, threads_, coresMax_);
        }
        ISAAC_THREAD_CERR << "Sorting " << std::distance(referenceSeedsBegin, referenceSeedsEnd) << " seeds done in " << (clock() - startSort) / 1000 << "ms" << std::endl;
        referenceSeedsBegin = referenceSeedsEnd;
    }

//    ISAAC_THREAD_CERR << "Removing Ns" << std::endl;
//    const size_t nCount = removeNs(seeds);
//    ISAAC_THREAD_CERR << "Removing Ns done (" << nCount << " removed - " << seeds.size() << " seeds remaining)" << std::endl;
}


/**
 * \brief Performs multiple passes If all seeds for unprocessed tiles don't fit in memory.
 *
 * \param unprocessedTiles [inout] All tiles are removed from the list upon return
 * \param tileClusterInfo  [in]   Cluster reads marked as complete are skipped.
 * \param foundMatches     [inout] Updated upon return.
 *
 */
template <typename KmerT>
void FindMatchesTransition::findMultiSeedMatches(
    const flowcell::Layout &flowcell,
    const std::vector<unsigned> &seedIndexList,
    flowcell::TileMetadataList &unprocessedTiles,
    alignment::matchFinder::TileClusterInfo &tileClusterInfo,
    SeedSource<KmerT> &seedSource,
    FoundMatchesMetadata &foundMatches)
{
    const unsigned seedLoaderOpenFileHandlesCount(inputLoadersMax_);

    alignment::SeedMetadataList seedMetadataList;
    BOOST_FOREACH(size_t index, seedIndexList)
    {
        seedMetadataList.push_back(flowcell.getSeedMetadataList().at(index));
    }
    ISAAC_THREAD_CERR << "Finding Multi-seed matches for " << seedMetadataList << " with repeat threshold: " <<
        repeatThreshold_ << std::endl;

    alignment::MatchFinder<KmerT> matchFinder(sortedReferenceMetadataList_, tempDirectory_,
                            unprocessedTiles, flowcell.getReadMetadataList(), flowcell.getSeedMetadataList(),
                            1,
                            ignoreNeighbors_, ignoreRepeats_,
                            repeatThreshold_, neighborhoodSizeThreshold_,
                            foundMatches.matchTally_, tileClusterInfo, threads_, coresMax_, tempSaversMax_,
                            standardOpenFileHandlesCount + seedLoaderOpenFileHandlesCount);

    flowcell::TileMetadataList currentTiles; currentTiles.reserve(unprocessedTiles.size());

    seedSource.initBuffers(unprocessedTiles, seedMetadataList);

    alignment::SeedMemoryManager<KmerT> seedMemoryManager(
        barcodeMetadataList_, flowcell.getReadMetadataList(), seedMetadataList, unprocessedTiles);

    while(!unprocessedTiles.empty())
    {
        if (!seedMemoryManager.selectTiles(unprocessedTiles, tileClusterInfo, matchFinder.getMaxTileCount(), tempSaversMax_, currentTiles))
        {
            BOOST_THROW_EXCEPTION(common::MemoryException("Insufficient memory to load seeds even for just one tile: " +
                boost::lexical_cast<std::string>(unprocessedTiles.back())));
        }

        {
            std::vector<alignment::Seed<KmerT> > seeds;
            seedMemoryManager.allocate(currentTiles, seeds);

            common::ScoopedMallocBlock  mallocBlock(memoryControl_);
            seedSource.generateSeeds(currentTiles, tileClusterInfo, seeds, mallocBlock);
            matchFinder.setTiles(currentTiles);

            ISAAC_THREAD_CERR << "Finding Exact multi-seed matches for " << seedMetadataList << " with repeat threshold: " <<
                repeatThreshold_ << std::endl;
            foundMatches.matchDistribution_.consolidate(
                matchFinder.findMatches(seeds.begin(), seedSource.getReferenceSeedBounds(), false, !neighborhoodSizeThreshold_));
            ISAAC_THREAD_CERR << "Finding Exact multi-seed matches done for " << seedMetadataList << std::endl;

            if (neighborhoodSizeThreshold_)
            {
                ISAAC_THREAD_CERR << "Finding Neighbor multi-seed matches for " << seedMetadataList << " with repeat threshold: " <<
                    repeatThreshold_ << std::endl;
                maskCompleteReadSeeds(
                    flowcell.getSeedMetadataList(), tileClusterInfo, seeds,
                    seedSource.getReferenceSeedBounds(), mallocBlock);
                foundMatches.matchDistribution_.consolidate(
                    matchFinder.findMatches(seeds.begin(), seedSource.getReferenceSeedBounds(), true, true));
                ISAAC_THREAD_CERR << "Finding Neighbor multi-seed matches done for " << seedMetadataList << std::endl;
            }
            std::vector<alignment::Seed<KmerT> >().swap(seeds);
        }
        currentTiles.clear();
    }

    ISAAC_THREAD_CERR << "Finding Multi-seed matches done for " << seedMetadataList << std::endl;

}

/**
 * \brief Finds matches for the lane. Updates foundMatches with match information and tile metadata identified during
 *        the processing.
 */
template <typename DataSourceT>
void FindMatchesTransition::findLaneMatches(
    const flowcell::Layout &flowcell,
    const unsigned lane,
    const flowcell::BarcodeMetadataList &laneBarcodes,
    flowcell::TileMetadataList &unprocessedTiles,
    DataSourceT &dataSource,
    demultiplexing::DemultiplexingStats &demultiplexingStats,
    FoundMatchesMetadata &foundMatches)
{
    if (!unprocessedTiles.empty())
    {
        alignment::matchFinder::TileClusterInfo tileClusterInfo(unprocessedTiles, clusterIdList_);
        ISAAC_THREAD_CERR << "Resolving barcodes for " << flowcell << " lane " << lane << std::endl;
        resolveBarcodes(
            flowcell, laneBarcodes, foundMatches.tileMetadataList_,
            dataSource,
            unprocessedTiles, tileClusterInfo, demultiplexingStats);
        ISAAC_THREAD_CERR << "Resolving barcodes done for " << flowcell << " lane " << lane << std::endl;

        const std::vector<std::vector<unsigned> > seedIndexListPerIteration = getSeedIndexListPerIteration(flowcell);

        while(!unprocessedTiles.empty())
        {
            flowcell::TileMetadataList thisPassTiles = findSingleSeedMatches(
                flowcell, seedIndexListPerIteration.at(0), 1 == seedIndexListPerIteration.size(),
                unprocessedTiles, tileClusterInfo, dataSource, demultiplexingStats,
                foundMatches);

            if (1 != seedIndexListPerIteration.size())
            {
                ISAAC_ASSERT_MSG(2 == seedIndexListPerIteration.size(), "2 sets of seeds expected");
                findMultiSeedMatches(flowcell, seedIndexListPerIteration.at(1), thisPassTiles, tileClusterInfo, dataSource, foundMatches);
                ISAAC_ASSERT_MSG(thisPassTiles.empty(), "Expected the findMultiSeedMatches to empty the list");
            }
        }
    }
}

inline bool orderByFlowcellLaneTile(const flowcell::TileMetadata& left, const flowcell::TileMetadata& right)
{
    return left.getFlowcellId() < right.getFlowcellId() ||
        (left.getFlowcellId() == right.getFlowcellId() &&
            (left.getLane() < right.getLane() ||
                (left.getLane() == right.getLane() && left.getTile() < right.getTile())));
}

inline const flowcell::TileMetadataList sortByFlowcellLaneTile(flowcell::TileMetadataList unsorted)
{
    std::sort(unsorted.begin(), unsorted.end(), orderByFlowcellLaneTile);
    return unsorted;
}

inline bool orderByFlowcellLane(const flowcell::BarcodeMetadata& left, const flowcell::BarcodeMetadata& right)
{
    return left.getFlowcellId() < right.getFlowcellId() ||
        (left.getFlowcellId() == right.getFlowcellId() &&
            (left.getLane() < right.getLane()));
}

inline const flowcell::BarcodeMetadataList sortByFlowcellLane(flowcell::BarcodeMetadataList unsorted)
{
    std::sort(unsorted.begin(), unsorted.end(), orderByFlowcellLane);
    return unsorted;
}

inline bool orderByFlowcell(const flowcell::Layout& left, const flowcell::Layout& right)
{
    return left.getFlowcellId() < right.getFlowcellId();
}

inline const flowcell::FlowcellLayoutList sortByFlowcell(flowcell::FlowcellLayoutList unsorted)
{
    std::sort(unsorted.begin(), unsorted.end(), orderByFlowcell);
    return unsorted;
}

/**
 * \brief Collect all barcodes belonging to the flowcell lane. Default barcode is placed at the beginning of the result list
 *
 * \return Returns subset of barcodeMetadataList or an empty list if none of the barcodes are mapped to a reference
 */
static flowcell::BarcodeMetadataList findFlowcellLaneBarcodes(
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const flowcell::Layout& flowcell,
    const unsigned lane)
{
    bool allBarcodesUnmapped = true;
    // put a placeholder for the unknown barcode in the beginning of the list
    flowcell::BarcodeMetadataList ret(1);
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        if (barcode.getFlowcellId() == flowcell.getFlowcellId() &&
            barcode.getLane() == lane)
        {
            ISAAC_THREAD_CERR << "adding " << barcode << std::endl;
            if (barcode.isDefault())
            {
                ISAAC_ASSERT_MSG(flowcell::BarcodeMetadata::INVALID_INDEX == ret.at(0).getIndex(), "More than one explicit specification for 'default' barcode within the group.");
                ret[0] = barcode;
            }
            else
            {
                ret.push_back(barcode);
            }
            allBarcodesUnmapped &= barcode.isUnmappedReference();
        }
    }

    ISAAC_ASSERT_MSG(flowcell::BarcodeMetadata::INVALID_INDEX != ret.at(0).getIndex(), "Missing default barcode specification");

    if (allBarcodesUnmapped)
    {
        ret.clear();
    }
    return ret;
}

template <typename DataSourceT>
void FindMatchesTransition::processFlowcellTiles(
    const flowcell::Layout& flowcell,
    DataSourceT &dataSource,
    demultiplexing::DemultiplexingStats &demultiplexingStats,
    FoundMatchesMetadata &foundMatches)
{
    TileSource &tileSource = dataSource;
    for (flowcell::TileMetadataList laneTiles = tileSource.discoverTiles(); !laneTiles.empty();
        laneTiles = tileSource.discoverTiles())
    {
        const unsigned lane = laneTiles[0].getLane();
        flowcell::BarcodeMetadataList laneBarcodes = findFlowcellLaneBarcodes(barcodeMetadataList_, flowcell, lane);
        if (laneBarcodes.empty())
        {
            ISAAC_THREAD_CERR << "Skipping flowcell " << flowcell.getFlowcellId() << " lane " << lane << " as none of the barcodes map to the reference" << std::endl;
        }
        else
        {
            BOOST_FOREACH(TileMetadata &tileMetadata, laneTiles)
            {
                foundMatches.addTile(tileMetadata);
                // this fixes the tile index to be correct in the context of the global tile list.
                tileMetadata = foundMatches.tileMetadataList_.back();
            }
            findLaneMatches(flowcell, lane, laneBarcodes, laneTiles, dataSource, demultiplexingStats, foundMatches);
        }
    }
}

template <typename KmerT>
void FindMatchesTransition::perform(FoundMatchesMetadata &foundMatches)
{
    FoundMatchesMetadata ret(tempDirectory_, barcodeMetadataList_, maxIterations_, sortedReferenceMetadataList_);
    demultiplexing::DemultiplexingStats demultiplexingStats(flowcellLayoutList_, barcodeMetadataList_);

    BOOST_FOREACH(const flowcell::Layout& flowcell, flowcellLayoutList_)
    {
        switch (flowcell.getFormat())
        {
            case flowcell::Layout::Bam:
            {
                BamSeedSource<KmerT> dataSource(
                    tempDirectory_,
                    availableMemory_,
                    clustersAtATimeMax_,
                    cleanupIntermediary_,
                    coresMax_, barcodeMetadataList_,
                    sortedReferenceMetadataList_, flowcell, threads_);
                processFlowcellTiles(flowcell, dataSource, demultiplexingStats, ret);
                break;
            }

            case flowcell::Layout::Fastq:
            {
                FastqSeedSource<KmerT> dataSource(
                    availableMemory_,
                    clustersAtATimeMax_,
                    allowVariableFastqLength_,
                    coresMax_, barcodeMetadataList_,
                    sortedReferenceMetadataList_, flowcell, threads_);
                processFlowcellTiles(flowcell, dataSource, demultiplexingStats, ret);
                break;
            }

            case flowcell::Layout::Bcl:
            {
                BclSeedSource<KmerT> dataSource(
                    ignoreMissingBcls_,
                    inputLoadersMax_, coresMax_, barcodeMetadataList_,
                    sortedReferenceMetadataList_, flowcell,
                    threads_);
                processFlowcellTiles(flowcell, dataSource, demultiplexingStats, ret);
                break;
            }
            case flowcell::Layout::BclBgzf:
            {
                BclBgzfSeedSource<KmerT> dataSource(
                    ignoreMissingBcls_,
                    inputLoadersMax_, coresMax_, barcodeMetadataList_,
                    sortedReferenceMetadataList_, flowcell,
                    threads_);
                processFlowcellTiles(flowcell, dataSource, demultiplexingStats, ret);
                break;
            }

            default:
            {
                ISAAC_ASSERT_MSG(false, "Unexpected flowcell format " << flowcell.getFormat());
                break;
            }
        }
    }

    dumpStats(demultiplexingStats, ret.tileMetadataList_);
    foundMatches.swap(ret);
}

void FindMatchesTransition::dumpStats(
    const demultiplexing::DemultiplexingStats &demultiplexingStats,
    const flowcell::TileMetadataList &tileMetadataList) const
{
    demultiplexing::DemultiplexingStatsXml statsXml;

    const unsigned maxLaneNumber = flowcell::getMaxLaneNumber(flowcellLayoutList_);
    for (unsigned lane = 1; lane <= maxLaneNumber; ++lane)
    {
        typedef std::map<std::string, demultiplexing::LaneBarcodeStats> SampleLaneBarcodeStats;
        typedef std::map<std::string, SampleLaneBarcodeStats> ProjectSampleLaneBarcodeStats;
        typedef std::map<std::string, ProjectSampleLaneBarcodeStats> FlowcellProjectSampleLaneBarcodeStats;
        FlowcellProjectSampleLaneBarcodeStats flowcellProjectSampleStats;
        BOOST_FOREACH(const flowcell::BarcodeMetadata& barcode, barcodeMetadataList_)
        {
            if (barcode.getLane() == lane)
            {
                // put one lane stat for each unknown barcode found.
                if (barcode.isUnknown())
                {
                    const flowcell::Layout& flowcell = flowcellLayoutList_.at(barcode.getFlowcellIndex());
                    statsXml.addFlowcellLane(flowcell, lane,
                                             demultiplexingStats.getLaneUnknwonBarcodeStat(barcode.getIndex()));
                }
                const demultiplexing::LaneBarcodeStats &stat = demultiplexingStats.getLaneBarcodeStat(barcode);
                flowcellProjectSampleStats[barcode.getFlowcellId()][barcode.getProject()][barcode.getSampleName()] += stat;
                flowcellProjectSampleStats[barcode.getFlowcellId()][barcode.getProject()]["all"] += stat;
                flowcellProjectSampleStats[barcode.getFlowcellId()]["all"]["all"] += stat;
                flowcellProjectSampleStats["all"][barcode.getProject()][barcode.getSampleName()] += stat;
                flowcellProjectSampleStats["all"][barcode.getProject()]["all"] += stat;
                flowcellProjectSampleStats["all"]["all"]["all"] += stat;
                statsXml.addLaneBarcode(barcode.getFlowcellId(), barcode.getProject(), barcode.getSampleName(), barcode.getName(), lane, stat);
            }
        }
        BOOST_FOREACH(const FlowcellProjectSampleLaneBarcodeStats::value_type &flowcellStats, flowcellProjectSampleStats)
        {
            BOOST_FOREACH(const ProjectSampleLaneBarcodeStats::value_type &projectStats, flowcellStats.second)
            {
                BOOST_FOREACH(const SampleLaneBarcodeStats::value_type &sampleStats, projectStats.second)
                {
                    statsXml.addLaneBarcode(flowcellStats.first, projectStats.first, sampleStats.first, "all", lane, sampleStats.second);
                }
            }
        }
    }

    std::ofstream os(demultiplexingStatsXmlPath_.string().c_str());
    if (!os) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: Unable to open file for writing: " + demultiplexingStatsXmlPath_.string()));
    }
    if (!(os << statsXml)) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: failed to store MatchFinder statistics in : " + demultiplexingStatsXmlPath_.string()));
    }
}


template void FindMatchesTransition::perform<oligo::ShortKmerType>(FoundMatchesMetadata &foundMatches);
template void FindMatchesTransition::perform<oligo::KmerType>(FoundMatchesMetadata &foundMatches);
template void FindMatchesTransition::perform<oligo::LongKmerType>(FoundMatchesMetadata &foundMatches);


} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
