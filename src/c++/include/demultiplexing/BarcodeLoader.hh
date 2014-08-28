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
 ** \file BarcodeLoader.hh
 **
 ** Helper class for loading barcode data from Bcl files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_BARCODE_LOADER_HH
#define iSAAC_DEMULTIPLEXING_BARCODE_LOADER_HH

#include <boost/iterator/counting_iterator.hpp>
#include <boost/mpl/equal_to.hpp>

#include "common/Debug.hh"
#include "common/Threads.hpp"

#include "demultiplexing/Barcode.hh"

#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/TileMetadata.hh"

#include "io/FileBufCache.hh"
#include "rta/BclMapper.hh"

namespace isaac
{
namespace demultiplexing
{


class BarcodeMemoryManager: boost::noncopyable
{
public:
    /// Determine how many tiles can have their barcoded loaded at the same time
    static bool selectTiles(
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

    static void allocate(const flowcell::TileMetadataList &tiles, std::vector<Barcode> &barcodes)
    {
        const unsigned long totalClusterCount = getTotalBarcodeCount(tiles);
        ISAAC_THREAD_CERR << "Allocating barcode storage for " << totalClusterCount << " barcodes" << std::endl;

        barcodes.clear();
        barcodes.resize(totalClusterCount);

        ISAAC_THREAD_CERR << "Allocating barcode storage done for " << totalClusterCount << " seeds" << std::endl;
    }

private:
    static bool seeIfFits(const flowcell::TileMetadataList &tiles)
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

    static unsigned getTotalBarcodeCount(const flowcell::TileMetadataList &tiles)
    {
        const unsigned long totalClusterCount = std::accumulate(tiles.begin(), tiles.end(), 0,
                                                             boost::bind(std::plus<unsigned>(), _1,
                                                                         boost::bind(&flowcell::TileMetadata::getClusterCount, _2)));
        return totalClusterCount;
    }
};


/**
 ** \brief Encapsulates the variables that are shared by all the threads while
 ** loading the barcodes.
 **/
template <typename ReaderT>
class ParallelBarcodeLoader : boost::noncopyable
{
public:
    /**
     ** \brief constructs an instance with all the required shorthands.
     **
     ** Note: all parameters are kept as references and it is the
     ** responsibility of the caller to ensure appropriate life time for the
     ** referenced variables.
     **/
    ParallelBarcodeLoader(
        const flowcell::TileMetadataList &tileMetadataList,
        const flowcell::Layout &flowcellLayout,
        const unsigned longestBclPath,
        std::vector<ReaderT> &threadReaders):
            flowcellLayout_(flowcellLayout)
            , maxClustersPerTile_(flowcell::getMaxTileClusters(tileMetadataList))
            , threadBclReaders_(threadReaders)
    {
        ISAAC_ASSERT_MSG(
            flowcell::Layout::Bcl == flowcellLayout_.getFormat() ||
            flowcell::Layout::BclBgzf == flowcellLayout_.getFormat(), "Unsupported flowcell format: " << flowcellLayout_.getFormat());

        ISAAC_ASSERT_MSG(MAX_BARCODE_LENGTH >= flowcellLayout_.getBarcodeLength(), "barcode cannot be longer than " << MAX_BARCODE_LENGTH << " bases");

        while(threadBclMappers_.size() < threadBclReaders_.size())
        {
            threadBclMappers_.push_back(
                new rta::SingleCycleBclMapper<ReaderT>(
                    maxClustersPerTile_,
                    longestBclPath, flowcell::Layout::Bcl != flowcellLayout_.getFormat(),
                    threadBclReaders_.at(threadBclMappers_.size())));
        }
    }


    void load(const unsigned unknownBarcodeIndex,
              std::vector<Barcode> &barcodes,
              std::vector<Barcode>::iterator &nextTileBarcodes,
              flowcell::TileMetadataList::const_iterator &nextTile,
              const flowcell::TileMetadataList::const_iterator tilesEnd,
              const unsigned threadNumber)
    {
        ISAAC_ASSERT_MSG(
            flowcell::Layout::Bcl == flowcellLayout_.getFormat() ||
            flowcell::Layout::BclBgzf == flowcellLayout_.getFormat(), "Only bcl barcode loading is supported");
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
                               boost::bind(&Barcode::constructFromTileBarcodeCluster, currentTile->getIndex(), unknownBarcodeIndex, _1));
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

private:
    // The mutex used to acquire the next tile and the destination of the seeds
    boost::mutex mutex_;
    const flowcell::Layout &flowcellLayout_;
    const unsigned maxClustersPerTile_;

    std::vector<ReaderT> &threadBclReaders_;
    boost::ptr_vector<rta::SingleCycleBclMapper<ReaderT> > threadBclMappers_;

    void loadTileCycle(
        rta::SingleCycleBclMapper<ReaderT> &bclMapper,
        std::vector<Barcode>::iterator destination,
        const flowcell::TileMetadata &tile,
        const unsigned cycle)
    {
        bclMapper.mapTileCycle(flowcellLayout_, tile, cycle);

        for (unsigned int clusterId = 0; tile.getClusterCount() > clusterId; ++clusterId)
        {
            char base = 0;
            bclMapper.get(clusterId, &base);
            Barcode &barcode = *destination++;

            const Kmer barcodeBase = (0 == base) ? 4 : (base & 3);
            barcode.setSequence((barcode.getSequence() << BITS_PER_BASE) | barcodeBase);
        }
    }
};

template <typename ReaderT>
class BarcodeLoader: boost::noncopyable
{
public:
    BarcodeLoader(
        common::ThreadVector &threads,
        const unsigned inputLoadersMax,
        const flowcell::TileMetadataList &allTilesMetadata,
        const flowcell::Layout &flowcellLayout,
        const unsigned longestBclPath,
        std::vector<ReaderT> &threadReaders):
            inputLoadersMax_(inputLoadersMax),
            threads_(threads),
            parallelBarcodeLoader_(allTilesMetadata, flowcellLayout, longestBclPath, threadReaders)
    {
    }

    /**
     * \brief resizes and fills barcodes with data.
     */
    void loadBarcodes(
        const unsigned unknownBarcodeIndex,
        const flowcell::TileMetadataList &tiles,
        std::vector<Barcode> &barcodes)
    {
        BarcodeMemoryManager::allocate(tiles, barcodes);
        // Start and execute the threads
        ISAAC_THREAD_CERR << "Loading data on " << inputLoadersMax_ << " threads" << std::endl;

        std::vector<flowcell::TileMetadata>::const_iterator nextTile = tiles.begin();
        std::vector<Barcode>::iterator nextTileBarcodes = barcodes.begin();
        threads_.execute(boost::bind(&ParallelBarcodeLoader<ReaderT>::load, &parallelBarcodeLoader_,
                                     unknownBarcodeIndex,
                                     boost::ref(barcodes),
                                     boost::ref(nextTileBarcodes),
                                     boost::ref(nextTile),
                                     tiles.end(), _1),
                         inputLoadersMax_);
    }

private:
    static unsigned getTotalBarcodeCount(const flowcell::TileMetadataList &tiles)
    {
        const unsigned long totalClusterCount = std::accumulate(tiles.begin(), tiles.end(), 0,
                                                             boost::bind(std::plus<unsigned>(), _1,
                                                                         boost::bind(&flowcell::TileMetadata::getClusterCount, _2)));
        return totalClusterCount;
    }

    const unsigned inputLoadersMax_;

    common::ThreadVector &threads_;
    ParallelBarcodeLoader<ReaderT> parallelBarcodeLoader_;
};

} // namespace demultiplexing
} // namespace isaac

#endif // #ifndef iSAAC_DEMULTIPLEXING_BARCODE_LOADER_HH
