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
 ** \file SeedMemoryManager.cpp
 **
 ** Deals with predicting and allocating memory for seeds.
 ** 
 ** \author Roman Petrovski
 **/
#include <algorithm>
#include <numeric>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include "alignment/SeedMemoryManager.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

template <typename KmerT>
SeedMemoryManager<KmerT>::SeedMemoryManager(
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const ReadMetadataList &readMetadataList,
    const SeedMetadataList &seedMetadataList,
    const flowcell::TileMetadataList &allTiles
    )
    : barcodeMetadataList_(barcodeMetadataList)
    , readMetadataList_(readMetadataList)
    , seedMetadataList_(seedMetadataList)
    , notFoundMatchesCount_()

{
    ISAAC_ASSERT_MSG(!readMetadataList_.empty(), "Empty readMetadataList is not allowed");
    ISAAC_ASSERT_MSG(!seedMetadataList_.empty(), "Empty seedMetadataList is not allowed");
}

template <typename KmerT>
bool SeedMemoryManager<KmerT>::seeIfFits(const TileMetadataList &tiles) const
{
    try
    {
        std::vector<SeedT> test;
        // * 2 is needed because the current implementation of parallelSort needs at least same size buffer for sorting
        // and a bit more...
        test.reserve(getTotalSeedCount(tiles) * 2 + 1024*1024*1024 / sizeof(SeedT));
        return true;
    }
    catch (std::bad_alloc &e)
    {
        // reset errno, to prevent misleading error messages when failing code does not set errno
        errno = 0;
    }
    return false;
}

template <typename KmerT>
bool SeedMemoryManager<KmerT>::selectTiles(TileMetadataList &unprocessedPool,
                             const matchFinder::TileClusterInfo &fragmentsToSkip,
                             const unsigned maxTilesAtATime,
                             const unsigned maxSavers,
                             TileMetadataList &selectedTiles)
{

    notFoundMatchesCount_ = getNotFoundMatchesCount(unprocessedPool, barcodeMetadataList_, readMetadataList_, fragmentsToSkip);
    selectedTiles.swap(unprocessedPool);
    {
        ISAAC_THREAD_CERR << "Determining the number of tiles that can be processed simultaneously..." << std::endl;
        while(!selectedTiles.empty() && (selectedTiles.size() > maxTilesAtATime || !seeIfFits(selectedTiles)))
        {
            // preserve the order. It is important.
            unprocessedPool.insert(unprocessedPool.begin(), selectedTiles.back());
            selectedTiles.pop_back();
        }
        if (selectedTiles.empty())
        {
            return false;
        }

        ISAAC_THREAD_CERR << "Determining the number of tiles that can be processed simultaneously done." << std::endl;

        if (!unprocessedPool.empty())
        {
            ISAAC_THREAD_CERR <<
                (selectedTiles.size() < maxTilesAtATime
                    ? "WARNING: will process tiles in parts due to the memory limit. "
                    : selectedTiles.size() == maxSavers
                      ? "WARNING: will process tiles in parts due to the parallel-save limit. "
                      : "WARNING: will process tiles in parts due to the open file handles limit. ") <<
                "This pass will process only " << selectedTiles.size() << " tiles" << std::endl;
        }
    }
    return true;
}

inline bool willLoadSeeds(
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const matchFinder::ClusterInfo &clusterInfo,
    const unsigned readIndex)
{
    return
        !clusterInfo.isReadComplete(readIndex) &&
        (!clusterInfo.isBarcodeSet() ||
            !barcodeMetadataList.at(clusterInfo.getBarcodeIndex()).isUnmappedReference());
}

template <typename KmerT>
const std::vector<std::vector<unsigned > > SeedMemoryManager<KmerT>::getNotFoundMatchesCount(
    const flowcell::TileMetadataList &unprocessedTiles,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const ReadMetadataList &readMetadataList,
    const matchFinder::TileClusterInfo &foundMatches) const
{
    ISAAC_ASSERT_MSG(unprocessedTiles.front().getIndex() <= unprocessedTiles.back().getIndex(),
                     "Expected tiles ordered by index");
    std::vector<std::vector<unsigned > > ret(readMetadataList.size(), std::vector<unsigned>(
        unprocessedTiles.back().getIndex() + 1));

    BOOST_FOREACH(const flowcell::ReadMetadata &readMetadata, readMetadataList_)
    {
        const unsigned readIndex = readMetadata.getIndex();
        BOOST_FOREACH(const flowcell::TileMetadata &tileMetadata, unprocessedTiles)
        {
            const unsigned tileIndex = tileMetadata.getIndex();
            const std::vector<matchFinder::ClusterInfo> &oneTileInfo = foundMatches.at(tileIndex);
            // match only clusters where no matches were found so far
            ISAAC_ASSERT_MSG(oneTileInfo.size() == tileMetadata.getClusterCount(), "allTiles and foundMatches geometries must match");
            ret.at(readIndex).at(tileIndex) =
                std::count_if(oneTileInfo.begin(), oneTileInfo.end(),
                              boost::bind(&willLoadSeeds,
                                          boost::ref(barcodeMetadataList), _1, readIndex));
        }
    }
    return ret;
}

template <typename KmerT>
unsigned long SeedMemoryManager<KmerT>::getTotalSeedCount(const TileMetadataList &tiles) const
{
    assert(!tiles.empty());

    unsigned long totalSeedCount = 0;
    BOOST_FOREACH(const flowcell::ReadMetadata &readMetadata, readMetadataList_)
    {
        const unsigned readIndex = readMetadata.getIndex();
        BOOST_FOREACH(const flowcell::TileMetadata &tileMetadata, tiles)
        {
            const unsigned tileIndex = tileMetadata.getIndex();
            // match only clusters where no matches were found so far
            const size_t notFoundCount = notFoundMatchesCount_.at(readIndex).at(tileIndex);
            BOOST_FOREACH(const SeedMetadata &seedMetadata, seedMetadataList_)
            {
                if (readIndex == seedMetadata.getReadIndex())
                {
                    totalSeedCount += notFoundCount * 2;
                }
            }
        }
    }
    return totalSeedCount;
}

template <typename KmerT>
void SeedMemoryManager<KmerT>::allocate(const TileMetadataList &tiles, std::vector<SeedT> &seeds) const
{
    // compute the total number of seeds that will be produced
    const unsigned long totalSeedCount = getTotalSeedCount(tiles);
    ISAAC_THREAD_CERR << "Allocating storage for "
                      << totalSeedCount
                      << " seeds (forward and reverse for "
                      << seedMetadataList_.size() << " seeds)" << std::endl;

    seeds.clear();
    seeds.resize(totalSeedCount);

    ISAAC_THREAD_CERR << "Allocating storage done for "
                      << totalSeedCount
                      << " seeds (forward and reverse for "
                      << seedMetadataList_.size() << " seeds)" << std::endl;
}

template class SeedMemoryManager<oligo::ShortKmerType>;
template class SeedMemoryManager<oligo::KmerType>;
template class SeedMemoryManager<oligo::LongKmerType>;

} // namespace alignment
} // namespace isaac
