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
 ** \file SeedGeneratorBase.cpp
 **
 ** Helper base class for seed generators
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

#include "alignment/SeedGeneratorBase.hh"
#include "common/Debug.hh"
#include "common/ParallelSort.hpp"
#include "common/SystemCompatibility.hh"

namespace isaac
{
namespace alignment
{

SeedGeneratorBase::SeedGeneratorBase(
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const flowcell::ReadMetadataList &readMetadataList,
    const std::vector<SeedMetadata> &seedMetadataList,
    const reference::SortedReferenceXmlList &sortedReferenceXmlList,
    const flowcell::TileMetadataList &tileMetadataList)
    : barcodeMetadataList_(barcodeMetadataList)
    , readMetadataList_(readMetadataList)
    , seedCounts_(getSeedCounts(readMetadataList_, seedMetadataList))
    , referenceTileReadFragmentCounts_(sortedReferenceXmlList.size(),
                                       std::vector<std::vector<unsigned> >(tileMetadataList.back().getIndex() + 1,
                                                                           std::vector<unsigned>(readMetadataList_.size())))
    , seedMetadataOrderedByFirstCycle_(orderSeedMetadataByFirstCycle(seedMetadataList))

    , nextTileSeedBegins_(sortedReferenceXmlList.size())
{
    ISAAC_ASSERT_MSG(!readMetadataList_.empty(), "Empty readMetadataList is not allowed");
    ISAAC_ASSERT_MSG(!seedMetadataList.empty(), "Empty seedMetadataList is not allowed");

    // referenceTileReadFragmentCounts_ needs to be sized so that getIndex stays within the vector.
    ISAAC_ASSERT_MSG(tileMetadataList.front().getIndex() <= tileMetadataList.back().getIndex(),
                     "Expected tiles ordered by index");
}

void SeedGeneratorBase::sortSeeds(
    std::vector<Seed> &seeds,
    common::ScoopedMallocBlock  &mallocBlock)
{
    // the Seeds buffer might turn out to be bigger than what we need if some clusters map to barcodes which
    // have unmapped references. This is not percieved to be the major scenario, so, some unused memory is acceptable
    ISAAC_ASSERT_MSG(seeds.end() >= getReferenceSeedBounds().back(), "Computed end is past the end of the reserved buffer");

    std::vector<Seed>::iterator referenceSeedsBegin = seeds.begin();
    BOOST_FOREACH(std::vector<Seed>::iterator referenceSeedsEnd, getReferenceSeedBounds())
    {
        ISAAC_THREAD_CERR << "Sorting " << referenceSeedsEnd - referenceSeedsBegin << " seeds" << std::endl;
        const clock_t startSort = clock();
        {
            common::ScoopedMallocBlockUnblock unblock(mallocBlock);
            // comparing the full kmer is required to push the N-seeds off to the very end.
            common::parallelSort(referenceSeedsBegin, referenceSeedsEnd, alignment::orderByKmerSeedIndex);
        }
        ISAAC_THREAD_CERR << "Sorting " << referenceSeedsEnd - referenceSeedsBegin << " seeds done in " << (clock() - startSort) / 1000 << "ms" << std::endl;
        referenceSeedsBegin = referenceSeedsEnd;
    }
}


inline bool firstCycleLess(const SeedMetadata &left, const SeedMetadata &right)
{
    return
        left.getReadIndex() < right.getReadIndex() ||
        (left.getReadIndex() == right.getReadIndex() && left.getOffset() < right.getOffset());
}

inline std::vector<SeedMetadata> SeedGeneratorBase::orderSeedMetadataByFirstCycle(
    std::vector<SeedMetadata> seedMetadataList)
{
    std::sort(seedMetadataList.begin(), seedMetadataList.end(), firstCycleLess);
    return seedMetadataList;
}

void SeedGeneratorBase::advanceToNextTile(const flowcell::TileMetadata &currentTile)
{
    // advance the marker to the destination of the seeds for the next tile
    const unsigned tileIndex = currentTile.getIndex();
    BOOST_FOREACH(const std::vector<std::vector<unsigned> > &reference, referenceTileReadFragmentCounts_)
    {
        const unsigned referenceIndex = &reference - &referenceTileReadFragmentCounts_.front();
        const std::vector<Seed>::iterator before = nextTileSeedBegins_.at(referenceIndex);
        BOOST_FOREACH(const flowcell::ReadMetadata &readMetadata, readMetadataList_)
        {
            const unsigned readIndex = readMetadata.getIndex();
            // forward and reverse seeds for this reference clusters excluding the 'complete ones' for this read
            nextTileSeedBegins_.at(referenceIndex) +=
                2 * reference.at(tileIndex).at(readIndex) * seedCounts_.at(readIndex);
        }
        ISAAC_THREAD_CERR << "RefIdx: " << referenceIndex << ", seeds: " << nextTileSeedBegins_.at(referenceIndex) - before << std::endl;
    }
}

std::vector<unsigned> SeedGeneratorBase::getSeedCounts(
    const std::vector<flowcell::ReadMetadata> &readMetadataList,
    const std::vector<SeedMetadata> &seedMetadataList) const
{
    std::vector<unsigned> seedCounts(readMetadataList.size(), 0);
    BOOST_FOREACH(const SeedMetadata &seedMetadata, seedMetadataList)
    {
        assert(seedCounts.size() > seedMetadata.getReadIndex());
        ++seedCounts[seedMetadata.getReadIndex()];
    }
    return seedCounts;
}

/**
 * \brief Computes the read fragment counts distribution for each tile having to be aligned to each reference
 * \param fragmentCounts [reference][tile][read] where the respective counts will be placed. The dimensions
 *                       must be pre-formatted
 */
void SeedGeneratorBase::getReferenceTileReadFragmentCounts(
    const flowcell::TileMetadataList &tiles,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const matchFinder::TileClusterInfo &tileClusterBarcode,
    FragmentCounts &fragmentCounts) const
{
    BOOST_FOREACH(std::vector<std::vector<unsigned> > &referenceFragmentCounts, fragmentCounts)
    {
        const unsigned referenceIndex = &referenceFragmentCounts - &fragmentCounts.front();

        BOOST_FOREACH(const flowcell::TileMetadata &tile, tiles)
        {
            std::vector<unsigned> &tileFragmentCounts = referenceFragmentCounts.at(tile.getIndex());
            BOOST_FOREACH(unsigned &readFragmentCount, tileFragmentCounts)
            {
                const unsigned readIndex = &readFragmentCount - &tileFragmentCounts.front();
                const std::vector<matchFinder::ClusterInfo> & thisTileReadsToSkip = tileClusterBarcode.at(tile.getIndex());

                readFragmentCount = 0;
                BOOST_FOREACH(const matchFinder::ClusterInfo &clusterInfo, thisTileReadsToSkip)
                {
                    ISAAC_ASSERT_MSG(clusterInfo.isBarcodeSet(), "Barcodes must be resolved at this point");
                    readFragmentCount +=
                        !clusterInfo.isReadComplete(readIndex) &&
                        barcodeMetadataList.at(clusterInfo.getBarcodeIndex()).getReferenceIndex() == referenceIndex;
                }
//                ISAAC_THREAD_CERR << "RefIdx: " << referenceIndex << ", readFragmentCount: " << readFragmentCount << std::endl;

            }
        }
    }
}

/**
 * \brief Format the vector of iterators so that each entry points at the first
 *        seed of each reference that will be loaded for the first tile
 */
void SeedGeneratorBase::reset(
    const flowcell::TileMetadataList &tiles, std::vector<Seed> &seeds,
    const matchFinder::TileClusterInfo &tileClusterBarcode)
{
    getReferenceTileReadFragmentCounts(tiles,
        barcodeMetadataList_, tileClusterBarcode, referenceTileReadFragmentCounts_);

    std::vector<Seed>::iterator referenceTilesBeginIterator = seeds.begin();
    nextTileSeedBegins_.clear();
    BOOST_FOREACH(const std::vector<std::vector<unsigned> >&reference, referenceTileReadFragmentCounts_)
    {
        nextTileSeedBegins_.push_back(referenceTilesBeginIterator);
        // use supplied tiles, not all the tiles in referenceTileReadFragmentCounts_!
        BOOST_FOREACH(const flowcell::TileMetadata &tile, tiles)
        {
            BOOST_FOREACH(const flowcell::ReadMetadata &read, readMetadataList_)
            {
                // reverse and forward seeds for each read that we decide to load.
                referenceTilesBeginIterator +=
                    2 * reference.at(tile.getIndex()).at(read.getIndex()) * seedCounts_.at(read.getIndex());
            }
        }
    }

}

} // namespace alignment
} // namespace isaac
