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


template <typename ReaderT, typename KmerT>
ParallelSeedLoader<ReaderT, KmerT>::ParallelSeedLoader(
    const bool ignoreMissingBcls,
    common::ThreadVector &threads,
    boost::ptr_vector<rta::SingleCycleBclMapper<ReaderT> > &threadBclMappers,
    const unsigned inputLoadersMax,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const flowcell::Layout &flowcellLayout,
    const std::vector<SeedMetadata> &seedMetadataList,
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
    const flowcell::TileMetadataList &tileMetadataList)
    : SeedGeneratorBase<KmerT>(barcodeMetadataList, flowcellLayout, seedMetadataList, sortedReferenceMetadataList, tileMetadataList)
    , inputLoadersMax_(inputLoadersMax)
    , seedCycles_(alignment::getAllSeedCycles(BaseT::flowcellLayout_.getReadMetadataList(), seedMetadataList))
    , threadDestinations_(inputLoadersMax, std::vector<typename std::vector<Seed<KmerT> >::iterator>(sortedReferenceMetadataList.size()))
    , threadCycleDestinations_(inputLoadersMax, std::vector<typename std::vector<Seed<KmerT> >::iterator>(sortedReferenceMetadataList.size()))
    , threads_(threads)
    , threadBclMappers_(threadBclMappers)
{
}

/**
 * \brief fills seeds with sorted ABCD permutation,
 *        N-containing seeds masked as poly-T with seed id ~0UL and moved to the back of each range.
 *
 * \param tiles               tiles to load
 * \param tileClusterBarcode  the index of the barcode for each cluster
 * \param seeds               storage for seeds
 */
template <typename ReaderT, typename KmerT>
void ParallelSeedLoader<ReaderT, KmerT>::loadSeeds(
    const flowcell::TileMetadataList &tiles,
    const matchFinder::TileClusterInfo &tileClusterBarcode,
    std::vector<Seed<KmerT> > &seeds,
    common::ScoopedMallocBlock  &mallocBlock)
{

    // Start and execute the threads
    ISAAC_THREAD_CERR << "Loading data on " << inputLoadersMax_ << " threads" << std::endl;

    std::vector<flowcell::TileMetadata>::const_iterator nextTile = tiles.begin();
    BaseT::reset(tiles, seeds, tileClusterBarcode);

    threads_.execute(boost::bind(&ParallelSeedLoader::load, this,
                                 boost::ref(tileClusterBarcode),
                                 boost::ref(nextTile),
                                 tiles.end(), _1),
                                 inputLoadersMax_);

    BaseT::sortSeeds(seeds, mallocBlock);
}

template <typename ReaderT, typename KmerT>
void ParallelSeedLoader<ReaderT, KmerT>::load(
    const matchFinder::TileClusterInfo &tileClusterBarcode,
    std::vector<flowcell::TileMetadata>::const_iterator &nextTile,
    const std::vector<flowcell::TileMetadata>::const_iterator tilesEnd,
    const unsigned threadNumber)
{
    // one iterator for each target range where seeds will be aligned against the same reference
    // using the class member vectors instead of stack to avoid vector memory allocation in the middle of seed loading
    std::vector<typename std::vector<Seed<KmerT> >::iterator> &thisThreadDestinationBegins = threadDestinations_.at(threadNumber);
    std::vector<typename std::vector<Seed<KmerT> >::iterator> &thisThreadCycleDestinationBegins = threadCycleDestinations_.at(threadNumber);
    boost::lock_guard<boost::mutex> lock(mutex_);
    while (tilesEnd != nextTile)
    {
        // acquire the next tile
        std::vector<flowcell::TileMetadata>::const_iterator currentTile = nextTile++;
        ISAAC_THREAD_CERR << "Loading tile seeds for " << *currentTile << std::endl;
        // acquire storage for forward and reverse kmer for all the seeds
        thisThreadDestinationBegins = BaseT::nextTileSeedBegins_;
        BaseT::advanceToNextTile(*currentTile);

        const unsigned tileIndex = currentTile->getIndex();
        {
            common::unlock_guard<boost::mutex> unlock(mutex_);
            std::vector<SeedMetadata>::const_iterator cycleSeedsBegin = BaseT::seedMetadataOrderedByFirstCycle_.begin();
            // cycles are guaranteed to belong to at least one of the seeds
            BOOST_FOREACH(const unsigned cycle, seedCycles_)
            {
                // find the first seed containing the cycle
                while (
                    cycleSeedsBegin->getOffset() +
                    BaseT::flowcellLayout_.getReadMetadataList().at(cycleSeedsBegin->getReadIndex()).getFirstCycle() +
                    cycleSeedsBegin->getLength() <= cycle)
                {
                    // advance each tile reference iterator to the first position of the first seed containing the cycle
                    BOOST_FOREACH(const std::vector<std::vector<unsigned> > &reference, BaseT::referenceTileReadFragmentCounts_)
                    {
                        const unsigned referenceIndex = &reference - &BaseT::referenceTileReadFragmentCounts_.front();

                        thisThreadDestinationBegins[referenceIndex] +=
                            2 * reference.at(tileIndex).at(cycleSeedsBegin->getReadIndex());
                    }
                    ++cycleSeedsBegin;
                }
                std::vector<SeedMetadata>::const_iterator cycleSeedsEnd = cycleSeedsBegin + 1;

                // find the first seed not containing the cycle
                while (
                    BaseT::seedMetadataOrderedByFirstCycle_.end() != cycleSeedsEnd &&
                    (cycleSeedsEnd->getOffset() +
                        BaseT::flowcellLayout_.getReadMetadataList().at(cycleSeedsEnd->getReadIndex()).getFirstCycle() <= cycle))
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

template <typename ReaderT, typename KmerT>
void ParallelSeedLoader<ReaderT, KmerT>::loadTileCycle(
    const matchFinder::TileClusterInfo &tileClusterBarcode,
    rta::SingleCycleBclMapper<ReaderT> &threadBclMapper,
    std::vector<typename std::vector<Seed<KmerT> >::iterator> &destinationBegins,
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

    threadBclMapper.mapTileCycle(BaseT::flowcellLayout_, tile, cycle);

    while (cycleSeedsEnd != cycleSeedsBegin)
    {
        for (unsigned int clusterId = 0; tile.getClusterCount() > clusterId; ++clusterId)
        {
            const unsigned barcodeIndex = clustersToDiscard.at(clusterId).getBarcodeIndex();
            const unsigned referenceIndex = BaseT::barcodeMetadataList_.at(barcodeIndex).getReferenceIndex();
            char base = 0;
            threadBclMapper.get(clusterId, &base);
            if (flowcell::BarcodeMetadata::UNMAPPED_REFERENCE_INDEX != referenceIndex)
            {
                if (!clustersToDiscard.at(clusterId).isReadComplete(readIndex))
                {
                    Seed<KmerT> & forwardSeed = *destinationBegins[referenceIndex]++;
                    Seed<KmerT> & reverseSeed = *destinationBegins[referenceIndex]++;
                    // skip those previously found to contain Ns
                    if (!forwardSeed.isNSeed())
                    {
                        if (!oligo::isBclN(base))
                        {
                            KmerT forward = forwardSeed.getKmer();
                            KmerT reverse = reverseSeed.getKmer();
                            const KmerT forwardBaseValue = base & oligo::BITS_PER_BASE_MASK;
                            const KmerT reverseBaseValue = (~forwardBaseValue) & oligo::BITS_PER_BASE_MASK;

                            forward <<= oligo::BITS_PER_BASE;
                            forward |= forwardBaseValue;
                            reverse >>= oligo::BITS_PER_BASE;
                            reverse |= (reverseBaseValue << (oligo::BITS_PER_BASE * oligo::KmerTraits<KmerT>::KMER_BASES - oligo::BITS_PER_BASE));

                            forwardSeed = Seed<KmerT>(forward, SeedId(tile.getIndex(), barcodeIndex, clusterId, cycleSeedsBegin->getIndex(), 0));
                            reverseSeed = Seed<KmerT>(reverse, SeedId(tile.getIndex(), barcodeIndex, clusterId, cycleSeedsBegin->getIndex(), 1));
                        }
                        else
                        {
                            // we can't have holes. The Ns must be stored in such a way that
                            // they will be easy to remove later (after sorting)
                            forwardSeed = makeNSeed<KmerT>(tile.getIndex(), barcodeIndex, clusterId, 0 == cycleSeedsBegin->getIndex());
                            reverseSeed = makeNSeed<KmerT>(tile.getIndex(), barcodeIndex, clusterId, 0 == cycleSeedsBegin->getIndex());
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

#include "rta/BclBgzfTileReader.hh"
template class isaac::alignment::ParallelSeedLoader<isaac::rta::BclBgzfTileReader, isaac::oligo::ShortKmerType>;
template class isaac::alignment::ParallelSeedLoader<isaac::rta::BclBgzfTileReader, isaac::oligo::KmerType>;
template class isaac::alignment::ParallelSeedLoader<isaac::rta::BclBgzfTileReader, isaac::oligo::LongKmerType>;

#include "rta/BclReader.hh"
template class isaac::alignment::ParallelSeedLoader<isaac::rta::BclReader, isaac::oligo::ShortKmerType>;
template class isaac::alignment::ParallelSeedLoader<isaac::rta::BclReader, isaac::oligo::KmerType>;
template class isaac::alignment::ParallelSeedLoader<isaac::rta::BclReader, isaac::oligo::LongKmerType>;
