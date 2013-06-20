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
 ** \file BarcodeLoader.cpp
 **
 ** Helper class for loading barcode data from Bcl files.
 **
 ** \author Roman Petrovski
 **/
#include <boost/iterator/counting_iterator.hpp>

#include "demultiplexing/BarcodeLoader.hh"

namespace isaac
{
namespace demultiplexing
{

BarcodeLoader::BarcodeLoader(
    const bool ignoreMissingBcls,
    common::ThreadVector &threads,
    const unsigned inputLoadersMax,
    const flowcell::TileMetadataList &allTilesMetadata,
    const flowcell::Layout &flowcellLayout,
    const flowcell::BarcodeMetadataList &barcodeMetadataList)
    : inputLoadersMax_(inputLoadersMax)
    , barcodeMetadataList_(barcodeMetadataList)
    , threads_(threads)
    , parallelBarcodeLoader_(ignoreMissingBcls, allTilesMetadata, flowcellLayout, barcodeMetadataList_, inputLoadersMax_)

{
}

static unsigned getTotalBarcodeCount(const flowcell::TileMetadataList &tiles)
{
    const unsigned long totalClusterCount = std::accumulate(tiles.begin(), tiles.end(), 0,
                                                         boost::bind(std::plus<unsigned>(), _1,
                                                                     boost::bind(&flowcell::TileMetadata::getClusterCount, _2)));
    return totalClusterCount;
}

void BarcodeLoader::allocate(const flowcell::TileMetadataList &tiles, std::vector<Barcode> &barcodes)
{
    const unsigned long totalClusterCount = getTotalBarcodeCount(tiles);
    ISAAC_THREAD_CERR << "Allocating barcode storage for " << totalClusterCount << " barcodes" << std::endl;

    barcodes.clear();
    barcodes.resize(totalClusterCount);

    ISAAC_THREAD_CERR << "Allocating barcode storage done for " << totalClusterCount << " seeds" << std::endl;
}


bool BarcodeLoader::seeIfFits(const flowcell::TileMetadataList &tiles)
{
    try
    {
        std::vector<Barcode> test;
        // * 2 is needed because the current implementation of parallelSort needs at least same size buffer for sorting
        // and a bit more...
        test.reserve(getTotalBarcodeCount(tiles) * 2 + 1024*1024*1024 / sizeof(Barcode));
        return true;
    }
    catch (std::bad_alloc &e)
    {
        // reset errno, to prevent misleading error messages when failing code does not set errno
        errno = 0;
    }
    return false;
}

bool BarcodeLoader::selectTiles(
    flowcell::TileMetadataList &unprocessedPool,
    flowcell::TileMetadataList &selectedTiles)
{
    selectedTiles.swap(unprocessedPool);
    {
        ISAAC_THREAD_CERR << "Barcode resolution: Determining the number of tiles that can be processed simultaneously..." << std::endl;
        while(!selectedTiles.empty() && !seeIfFits(selectedTiles))
        {
            unprocessedPool.push_back(selectedTiles.back());
            selectedTiles.pop_back();
        }
        if (selectedTiles.empty())
        {
            return false;
        }

        ISAAC_THREAD_CERR << "Barcode resolution: Determining the number of tiles that can be processed simultaneously done." << std::endl;

        if (!unprocessedPool.empty())
        {
            ISAAC_THREAD_CERR << "WARNING: will resolve barcodes in parts due to the memory limit. "
                "This pass will process only " << selectedTiles.size() << " tiles" << std::endl;
        }
    }
    return true;
}


/**
 * \brief Loads barcodes
 */
void BarcodeLoader::loadBarcodes(const flowcell::TileMetadataList &tiles, std::vector<Barcode> &barcodes)
{
    // Start and execute the threads
    ISAAC_THREAD_CERR << "Loading data on " << inputLoadersMax_ << " threads" << std::endl;

    std::vector<flowcell::TileMetadata>::const_iterator nextTile = tiles.begin();
    std::vector<Barcode>::iterator nextTileBarcodes = barcodes.begin();
    threads_.execute(boost::bind(&ParallelBarcodeLoader::load, &parallelBarcodeLoader_,
                                 boost::ref(barcodes),
                                 boost::ref(nextTileBarcodes),
                                 boost::ref(nextTile),
                                 tiles.end(), _1),
                     inputLoadersMax_);
}

ParallelBarcodeLoader::ParallelBarcodeLoader(
    const bool ignoreMissingBcls,
    const flowcell::TileMetadataList &tileMetadataList,
    const flowcell::Layout &flowcellLayout,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const unsigned maxThreads)
    : tileMetadataList_(tileMetadataList)
    , flowcellLayout_(flowcellLayout)
    , unknownBarcodeIndex_(barcodeMetadataList.at(0).getIndex())
{
    if (flowcell::Layout::Bcl == flowcellLayout_.getFormat() ||
        flowcell::Layout::BclGz == flowcellLayout_.getFormat())
    {
        ISAAC_ASSERT_MSG(MAX_BARCODE_LENGTH >= flowcellLayout_.getBarcodeLength(), "barcode cannot be longer than " << MAX_BARCODE_LENGTH << " bases");
        ISAAC_ASSERT_MSG(barcodeMetadataList.size(), "Barcode list must be not empty");
        ISAAC_ASSERT_MSG(barcodeMetadataList.at(0).isDefault(), "The very first barcode must be the 'unknown indexes or no index' one");

        const unsigned highestTileNumber = std::max_element(tileMetadataList_.begin(), tileMetadataList_.end(),
                                                            boost::bind(&flowcell::TileMetadata::getTile, _1)<
                                                                boost::bind(&flowcell::TileMetadata::getTile, _2))->getTile();
        const unsigned insanelyHighCycleNumber = 9999;

        boost::filesystem::path longestBclFilePath;
        flowcellLayout_.getBclFilePath(highestTileNumber, 1, insanelyHighCycleNumber, longestBclFilePath);

        const unsigned maxClusterCount = std::max_element(tileMetadataList_.begin(), tileMetadataList_.end(),
                                                          boost::bind(&flowcell::TileMetadata::getClusterCount, _1)<
                                                          boost::bind(&flowcell::TileMetadata::getClusterCount, _2))->getClusterCount();

        while(threadBclMappers_.size() < maxThreads)
        {
            threadBclMappers_.push_back(new io::SingleCycleBclMapper(ignoreMissingBcls, maxClusterCount));
            threadBclMappers_.back().reserveBuffers(longestBclFilePath.string().size(), flowcell::Layout::BclGz == flowcellLayout.getFormat());
        }
    }
}

void ParallelBarcodeLoader::load(std::vector<Barcode> &barcodes,
                              std::vector<Barcode>::iterator &nextTileBarcodes,
                              flowcell::TileMetadataList::const_iterator &nextTile,
                              const flowcell::TileMetadataList::const_iterator tilesEnd,
                              const unsigned threadNumber)
{
    ISAAC_ASSERT_MSG(
        flowcell::Layout::Bcl == flowcellLayout_.getFormat() ||
        flowcell::Layout::BclGz == flowcellLayout_.getFormat(), "Only bcl barcode loading is supported");
    boost::lock_guard<boost::mutex> lock(mutex_);
    while (tilesEnd != nextTile)
    {
        flowcell::TileMetadataList::const_iterator currentTile = nextTile++;
        std::vector<Barcode>::iterator destinationBegin = nextTileBarcodes;
        nextTileBarcodes += currentTile->getClusterCount();

        const std::vector<Barcode>::const_iterator  destinationEnd = nextTileBarcodes;
        ISAAC_ASSERT_MSG(destinationEnd <= barcodes.end(), "Computed end is past the end of the reserved buffer");

        {
            common::unlock_guard<boost::mutex> unlock(mutex_);

            ISAAC_THREAD_CERR << "Formatting tile barcodes for " << *currentTile << std::endl;
            std::transform(boost::counting_iterator<unsigned>(0),
                           boost::counting_iterator<unsigned>(currentTile->getClusterCount()),
                           destinationBegin,
                           boost::bind(&Barcode::constructFromTileBarcodeCluster, currentTile->getIndex(), unknownBarcodeIndex_, _1));
            ISAAC_THREAD_CERR << "Formatting tile barcodes done for " << *currentTile << std::endl;

            ISAAC_THREAD_CERR << "Loading tile barcodes for " << *currentTile << std::endl;
            BOOST_FOREACH(const unsigned cycle, flowcellLayout_.getBarcodeCycles())
            {
                loadTileCycle(threadBclMappers_.at(threadNumber), destinationBegin, *currentTile, cycle);
            }
            ISAAC_THREAD_CERR << "Loading tile barcodes done for " << *currentTile << std::endl;
        }
    }
}

void ParallelBarcodeLoader::loadTileCycle(
    io::SingleCycleBclMapper &threadBclMapper,
    std::vector<Barcode>::iterator destination,
    const flowcell::TileMetadata &tile,
    const unsigned cycle)
{
    threadBclMapper.mapTileCycle(flowcellLayout_, tile, cycle);

    for (unsigned int clusterId = 0; tile.getClusterCount() > clusterId; ++clusterId)
    {
        char base = 0;
        threadBclMapper.get(clusterId, &base);
        Barcode &barcode = *destination++;

        const Kmer barcodeBase = (0 == base) ? 4 : (base & 3);
        barcode.setSequence((barcode.getSequence() << BITS_PER_BASE) | barcodeBase);
    }
}

} // namespace demultiplexing
} // namespace isaac

