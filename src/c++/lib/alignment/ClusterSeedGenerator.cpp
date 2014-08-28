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
 ** \file ClusterSeedGenerator.cpp
 **
 ** Component to generate the seeds from a block of sequentially-stored bcl clusters
 ** 
 ** \author Roman Petrovski
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

#include "alignment/BclClusters.hh"
#include "alignment/ClusterSeedGenerator.hh"
#include "common/SystemCompatibility.hh"

namespace isaac
{
namespace alignment
{

template <typename KmerT>
ClusterSeedGenerator<KmerT>::ClusterSeedGenerator(
    common::ThreadVector &threads,
    const unsigned computeThreadsMax,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const flowcell::Layout &flowcellLayout,
    const std::vector<SeedMetadata> &seedMetadataList,
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
    const flowcell::TileMetadataList &willRequestTiles,
    const BclClusters &clusters,
    const flowcell::TileMetadataList &loadedTiles)
    : SeedGeneratorBase<KmerT>(barcodeMetadataList, flowcellLayout, seedMetadataList, sortedReferenceMetadataList, willRequestTiles)
    , clusters_(clusters)
    , loadedTiles_(loadedTiles)
    , computeThreadsMax_(computeThreadsMax)
    , threadDestinations_(computeThreadsMax_, std::vector<typename std::vector<Seed<KmerT> >::iterator>(sortedReferenceMetadataList.size()))
    , threads_(threads)
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
template <typename KmerT>
void ClusterSeedGenerator<KmerT>::generateSeeds(
    const flowcell::TileMetadataList &tiles,
    const matchFinder::TileClusterInfo &tileClusterBarcode,
    std::vector<Seed<KmerT> > &seeds,
    common::ScoopedMallocBlock  &mallocBlock)
{

    // Start and execute the threads
    ISAAC_THREAD_CERR << "Loading data on " << computeThreadsMax_ << " threads" << std::endl;

    std::vector<flowcell::TileMetadata>::const_iterator nextTile = tiles.begin();
    BaseT::reset(tiles, seeds, tileClusterBarcode);

    threads_.execute(boost::bind(&ClusterSeedGenerator::generateThread, this,
                                 boost::ref(tileClusterBarcode),
                                 boost::ref(nextTile),
                                 tiles.begin(), tiles.end(), _1),
                                 computeThreadsMax_);

    BaseT::sortSeeds(seeds, mallocBlock);
}

template <typename KmerT>
void ClusterSeedGenerator<KmerT>::generateThread(
    const matchFinder::TileClusterInfo &tileClusterBarcode,
    std::vector<flowcell::TileMetadata>::const_iterator &nextTile,
    const std::vector<flowcell::TileMetadata>::const_iterator tilesBegin,
    const std::vector<flowcell::TileMetadata>::const_iterator tilesEnd,
    const unsigned threadNumber)
{
    // one iterator for each target range where seeds will be aligned against the same reference
    // using the class member vectors instead of stack to avoid vector memory allocation in the middle of seed loading
    std::vector<typename std::vector<Seed<KmerT> >::iterator> &thisThreadDestinationBegins = threadDestinations_.at(threadNumber);
    boost::lock_guard<boost::mutex> lock(mutex_);
    while (tilesEnd != nextTile)
    {
        // acquire the next tile
        flowcell::TileMetadataList::const_iterator currentTile = nextTile++;
        ISAAC_THREAD_CERR << "Generating seeds for " << *currentTile << std::endl;
        // acquire storage for forward and reverse kmer for all the seeds
        thisThreadDestinationBegins = BaseT::nextTileSeedBegins_;

        BaseT::advanceToNextTile(*currentTile);

        const unsigned tileIndex = currentTile->getIndex();
        const std::vector<matchFinder::ClusterInfo> &clustersToDiscard = tileClusterBarcode.at(tileIndex);
        ISAAC_ASSERT_MSG(currentTile->getClusterCount() == clustersToDiscard.size(), "Found matches from a wrong tile/read");
        {
            common::unlock_guard<boost::mutex> unlock(mutex_);

            unsigned firstClusterOffset = 0;
            // sum all the tile clusters before the one we're processing ot find the offset of our first cluster in the entire loaded clusters thing.
            for (flowcell::TileMetadataList::const_iterator preceedingTileIt = loadedTiles_.begin();
                ; ++preceedingTileIt)
            {
                ISAAC_ASSERT_MSG(loadedTiles_.end() != preceedingTileIt, "Loaded tiles list does not contain the tile we're about to generate seeds for");
                ISAAC_ASSERT_MSG(currentTile->getFlowcellIndex() == preceedingTileIt->getFlowcellIndex(),
                                 "Expecting all tiles to belong to the same flowcell");
                if (currentTile->getLane() == preceedingTileIt->getLane() &&
                    currentTile->getTile() == preceedingTileIt->getTile())
                {
                    break;
                }
                firstClusterOffset += preceedingTileIt->getClusterCount();
            }

            BclClusters::const_iterator clusterIt = clusters_.cluster(firstClusterOffset);
            BclClusters::const_iterator clustersEnd = clusters_.cluster(firstClusterOffset + currentTile->getClusterCount());
            unsigned long long totalSeeds = 0;
            unsigned long long nSeeds = 0;
            for (unsigned clusterId = 0; clustersEnd != clusterIt; ++clusterId, clusterIt += clusters_.getClusterLength())
            {
                const unsigned barcodeIndex = clustersToDiscard.at(clusterId).getBarcodeIndex();
                const unsigned referenceIndex = BaseT::barcodeMetadataList_.at(barcodeIndex).getReferenceIndex();

                if (flowcell::BarcodeMetadata::UNMAPPED_REFERENCE_INDEX != referenceIndex)
                {
                    BOOST_FOREACH(const SeedMetadata &seedMetadata, BaseT::seedMetadataOrderedByFirstCycle_)
                    {
                        const unsigned readIndex = seedMetadata.getReadIndex();
                        if (!clustersToDiscard.at(clusterId).isReadComplete(readIndex))
                        {
                            Seed<KmerT> & forwardSeed = *thisThreadDestinationBegins[referenceIndex]++;
                            Seed<KmerT> & reverseSeed = *thisThreadDestinationBegins[referenceIndex]++;

                            forwardSeed = Seed<KmerT>(0, SeedId(currentTile->getIndex(), barcodeIndex, clusterId, seedMetadata.getIndex(), 0));
                            reverseSeed = Seed<KmerT>(0, SeedId(currentTile->getIndex(), barcodeIndex, clusterId, seedMetadata.getIndex(), 1));

                            unsigned len = oligo::KmerTraits<KmerT>::KMER_BASES;
                            for (BclClusters::const_iterator baseIt =
                                    clusterIt + seedMetadata.getOffset() + BaseT::flowcellLayout_.getReadMetadataList().at(readIndex).getOffset();
                                len; --len, ++baseIt)
                            {
                                const unsigned char base = *baseIt;
                                if (!oligo::isBclN(base))
                                {
                                    const KmerT forwardBaseValue = base & oligo::BITS_PER_BASE_MASK;
                                    const KmerT reverseBaseValue = (~forwardBaseValue) & oligo::BITS_PER_BASE_MASK;

                                    forwardSeed.kmer() <<= oligo::BITS_PER_BASE;
                                    forwardSeed.kmer() |= forwardBaseValue;
                                    reverseSeed.kmer() >>= oligo::BITS_PER_BASE;
                                    reverseSeed.kmer() |= (reverseBaseValue << (oligo::BITS_PER_BASE * oligo::KmerTraits<KmerT>::KMER_BASES - oligo::BITS_PER_BASE));

                                    if (len == 1)
                                    {
                                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, "fw:" << forwardSeed);
                                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, "rv:" << reverseSeed);
                                    }
                                }
                                else
                                {
                                    // we can't have holes. The Ns must be stored in such a way that
                                    // they will be easy to remove later (after sorting)
                                    forwardSeed = makeNSeed<KmerT>(currentTile->getIndex(), barcodeIndex, clusterId, 0 == seedMetadata.getIndex());
                                    reverseSeed = makeNSeed<KmerT>(currentTile->getIndex(), barcodeIndex, clusterId, 0 == seedMetadata.getIndex());
                                    ++nSeeds;
                                    break;
                                }
                            }
                            ++totalSeeds;
                        }
                    }
                }
            }
            ISAAC_THREAD_CERR << "Generating seeds done for " << *currentTile << " generated " <<
                (totalSeeds - nSeeds) * 2 << " normal and " << nSeeds * 2 << " n-seeds" << std::endl;
        }
    }
}

template class ClusterSeedGenerator<oligo::ShortKmerType>;
template class ClusterSeedGenerator<oligo::KmerType>;
template class ClusterSeedGenerator<oligo::LongKmerType>;

} // namespace alignment
} // namespace isaac
