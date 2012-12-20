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
 ** \file SeedLoader.cpp
 **
 ** Component to load the seeds in parallel
 ** 
 ** \author Come Raczy
 **/
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cstring>
#include <cerrno>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/thread.hpp>
#include <boost/foreach.hpp>

#include "alignment/SeedLoader.hh"
#include "common/Debug.hh"
#include "common/ParallelSort.hpp"
#include "common/SystemCompatibility.hh"

namespace isaac
{
namespace alignment
{


ParallelSeedLoader::ParallelSeedLoader(
    const bool ignoreMissingBcls,
    common::ThreadVector &threads,
    const unsigned inputLoadersMax,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const flowcell::ReadMetadataList &readMetadataList,
    const std::vector<SeedMetadata> &seedMetadataList,
    const reference::SortedReferenceXmlList &sortedReferenceXmlList,
    const flowcell::TileMetadataList &tileMetadataList)
    : SeedGeneratorBase(barcodeMetadataList, readMetadataList, seedMetadataList, sortedReferenceXmlList, tileMetadataList)
    , inputLoadersMax_(inputLoadersMax)
    , seedCycles_(discoverSeedCycles())
    , threadDestinations_(inputLoadersMax, std::vector<std::vector<Seed>::iterator>(sortedReferenceXmlList.size()))
    , threadCycleDestinations_(inputLoadersMax, std::vector<std::vector<Seed>::iterator>(sortedReferenceXmlList.size()))
    , threads_(threads)
{
    const boost::filesystem::path longestBaseCallsPath =  flowcell::getLongestBaseCallsPath(tileMetadataList);
    const unsigned highestTileNumber = std::max_element(tileMetadataList.begin(), tileMetadataList.end(),
                                                        boost::bind(&flowcell::TileMetadata::getTile, _1)<
                                                        boost::bind(&flowcell::TileMetadata::getTile, _2))->getTile();
    boost::filesystem::path longestBclFilePath;
    flowcell::Layout::getBclFilePath(highestTileNumber, 1, longestBaseCallsPath, seedCycles_.back(),
                                     flowcell::TileMetadata::GzCompression, longestBclFilePath);

    const unsigned maxClusterCount = std::max_element(tileMetadataList.begin(), tileMetadataList.end(),
                                                       boost::bind(&flowcell::TileMetadata::getClusterCount, _1)<
                                                       boost::bind(&flowcell::TileMetadata::getClusterCount, _2))->getClusterCount();
    const bool compressedBclsFound =
        tileMetadataList.end() != std::find_if(tileMetadataList.begin(), tileMetadataList.end(),
                                                boost::bind(&flowcell::TileMetadata::getCompression, _1) !=
                                                    flowcell::TileMetadata::NoCompression);
    while(threadBclMappers_.size() < inputLoadersMax_)
    {
        threadBclMappers_.push_back(new io::SingleCycleBclMapper(ignoreMissingBcls, maxClusterCount));
        threadBclMappers_.back().reserveBuffers(longestBclFilePath.string().size(), compressedBclsFound);
    }
}

/**
 * \brief fills seeds with sorted ABCD permuation,
 *        N-containing seeds masked as poly-T with seed id ~0UL and moved to the back of each range.
 *
 * \param tiles               tiles to load
 * \param tileClusterBarcode  the index of the barcode for each cluster
 * \param seeds               storage for seeds
 */
void ParallelSeedLoader::loadSeeds(
    const flowcell::TileMetadataList &tiles,
    const matchFinder::TileClusterInfo &tileClusterBarcode,
    std::vector<Seed> &seeds,
    common::ScoopedMallocBlock  &mallocBlock)
{

    // Start and execute the threads
    ISAAC_THREAD_CERR << "Loading data on " << inputLoadersMax_ << " threads" << std::endl;

    std::vector<flowcell::TileMetadata>::const_iterator nextTile = tiles.begin();
    reset(tiles, seeds, tileClusterBarcode);

    threads_.execute(boost::bind(&ParallelSeedLoader::load, this,
                                 boost::ref(tileClusterBarcode),
                                 boost::ref(nextTile),
                                 tiles.end(), _1),
                                 inputLoadersMax_);

    sortSeeds(seeds, mallocBlock);
}

std::vector<unsigned> ParallelSeedLoader::discoverSeedCycles() const
{
    std::vector<unsigned> ret;
    BOOST_FOREACH(const SeedMetadata &seedMetadata, seedMetadataOrderedByFirstCycle_)
    {
        const unsigned seedFirstCycle =
            seedMetadata.getOffset() + readMetadataList_.at(seedMetadata.getReadIndex()).getFirstCycle();
        for(unsigned cycle = seedFirstCycle; cycle < seedFirstCycle + seedMetadata.getLength(); ++cycle)
        {
            if (ret.empty() || ret.back() < cycle)
            {
                ret.push_back(cycle);
            }
        }
    }
    return ret;
}

void ParallelSeedLoader::load(
    const matchFinder::TileClusterInfo &tileClusterBarcode,
    std::vector<flowcell::TileMetadata>::const_iterator &nextTile,
    const std::vector<flowcell::TileMetadata>::const_iterator tilesEnd,
    const unsigned threadNumber)
{
    // one iterator for each target range where seeds will be aligned against the same reference
    // using the class member vectors instead of stack to avoid vector memory allocation in the middle of seed loading
    std::vector<std::vector<Seed>::iterator> &thisThreadDestinationBegins = threadDestinations_.at(threadNumber);
    std::vector<std::vector<Seed>::iterator> &thisThreadCycleDestinationBegins = threadCycleDestinations_.at(threadNumber);
    boost::lock_guard<boost::mutex> lock(mutex_);
    while (tilesEnd != nextTile)
    {
        // acquire the next tile
        std::vector<flowcell::TileMetadata>::const_iterator currentTile = nextTile++;
        ISAAC_THREAD_CERR << "Loading tile seeds for " << *currentTile << std::endl;
        // acquire storage for forward and reverse kmer for all the seeds
        thisThreadDestinationBegins = nextTileSeedBegins_;
        advanceToNextTile(*currentTile);

        const unsigned tileIndex = currentTile->getIndex();
        {
            common::unlock_guard<boost::mutex> unlock(mutex_);
            std::vector<SeedMetadata>::const_iterator cycleSeedsBegin = seedMetadataOrderedByFirstCycle_.begin();
            // cycles are guaranteed to belong to at least one of the seeds
            BOOST_FOREACH(const unsigned cycle, seedCycles_)
            {
                // find the first seed containing the cycle
                while (
                    cycleSeedsBegin->getOffset() +
                    readMetadataList_.at(cycleSeedsBegin->getReadIndex()).getFirstCycle() +
                    cycleSeedsBegin->getLength() <= cycle)
                {
                    // advance each tile reference iterator to the first position of the first seed containing the cycle
                    BOOST_FOREACH(const std::vector<std::vector<unsigned> > &reference, referenceTileReadFragmentCounts_)
                    {
                        const unsigned referenceIndex = &reference - &referenceTileReadFragmentCounts_.front();

                        thisThreadDestinationBegins[referenceIndex] +=
                            2 * reference.at(tileIndex).at(cycleSeedsBegin->getReadIndex());
                    }
                    ++cycleSeedsBegin;
                }
                std::vector<SeedMetadata>::const_iterator cycleSeedsEnd = cycleSeedsBegin + 1;

                // find the first seed not containing the cycle
                while (
                    seedMetadataOrderedByFirstCycle_.end() != cycleSeedsEnd &&
                    (cycleSeedsEnd->getOffset() +
                        readMetadataList_.at(cycleSeedsEnd->getReadIndex()).getFirstCycle() <= cycle))
                {
                    ++cycleSeedsEnd;
                }

                // this call messes up thisThreadCycleDestinationBegins so, save it.
                thisThreadCycleDestinationBegins = thisThreadDestinationBegins;
                loadTileCycle(tileClusterBarcode, threadBclMappers_.at(threadNumber),
                              thisThreadCycleDestinationBegins, *currentTile, cycle, cycleSeedsBegin, cycleSeedsEnd);
            }
            ISAAC_THREAD_CERR << "Loading tile seeds done for " << *currentTile << std::endl;
        }
    }
}

void ParallelSeedLoader::loadTileCycle(
    const matchFinder::TileClusterInfo &tileClusterBarcode,
    io::SingleCycleBclMapper &threadBclMapper,
    std::vector<std::vector<Seed>::iterator> &destinationBegins,
    const flowcell::TileMetadata &tile,
    const unsigned cycle,
    std::vector<SeedMetadata>::const_iterator cycleSeedsBegin,
    const std::vector<SeedMetadata>::const_iterator cycleSeedsEnd)
{
    ISAAC_ASSERT_MSG(cycleSeedsEnd > cycleSeedsBegin, "Seed list cannot be empty");
    ISAAC_ASSERT_MSG((cycleSeedsEnd -1)->getReadIndex() == cycleSeedsBegin->getReadIndex(), "All seeds must belong to the same read");
    const unsigned readIndex = cycleSeedsBegin->getReadIndex();

    const std::vector<matchFinder::ClusterInfo> &clustersToDiscard = tileClusterBarcode.at(tile.getIndex());
    ISAAC_ASSERT_MSG(tile.getClusterCount() == clustersToDiscard.size(), "Found matches from a wrong tile/read");

    threadBclMapper.mapTileCycle(tile, cycle);

    while (cycleSeedsEnd != cycleSeedsBegin)
    {
        for (unsigned int clusterId = 0; tile.getClusterCount() > clusterId; ++clusterId)
        {
            const unsigned barcodeIndex = clustersToDiscard.at(clusterId).getBarcodeIndex();
            const unsigned referenceIndex = barcodeMetadataList_.at(barcodeIndex).getReferenceIndex();
            char base = 0;
            threadBclMapper.get(clusterId, &base);
            if (flowcell::BarcodeMetadata::UNMAPPED_REFERENCE_INDEX != referenceIndex)
            {
                if (!clustersToDiscard.at(clusterId).isReadComplete(readIndex))
                {
                    Seed & forwardSeed = *destinationBegins[referenceIndex]++;
                    Seed & reverseSeed = *destinationBegins[referenceIndex]++;
                    // skip those previously found to contain Ns
                    if (!forwardSeed.isNSeed())
                    {
                        if (!oligo::isBclN(base))
                        {
                            oligo::Kmer forward = forwardSeed.getKmer();
                            oligo::Kmer reverse = reverseSeed.getKmer();
                            const oligo::Kmer forwardBaseValue = base & 3;
                            const oligo::Kmer reverseBaseValue = (~forwardBaseValue) & 3;

                            forward <<= 2;
                            forward |= forwardBaseValue;
                            reverse >>= 2;
                            reverse |= (reverseBaseValue << (2 * oligo::kmerLength - 2));

                            forwardSeed = Seed(forward, SeedId(tile.getIndex(), barcodeIndex, clusterId, cycleSeedsBegin->getIndex(), 0));
                            reverseSeed = Seed(reverse, SeedId(tile.getIndex(), barcodeIndex, clusterId, cycleSeedsBegin->getIndex(), 1));
                        }
                        else
                        {
                            // we can't have holes. The Ns must be stored in such a way that
                            // they will be easy to remove later (after sorting)
                            forwardSeed = makeNSeed(tile.getIndex(), barcodeIndex, clusterId, 0 == cycleSeedsBegin->getIndex());
                            reverseSeed = makeNSeed(tile.getIndex(), barcodeIndex, clusterId, 0 == cycleSeedsBegin->getIndex());
                        }
                    }
                }
            }
        }
        ++cycleSeedsBegin;
    }
}


} // namespace alignment
} // namespace isaac
