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
 ** \file BarcodeLoader.hh
 **
 ** Helper class for loading barcode data from Bcl files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_BARCODE_LOADER_HH
#define iSAAC_DEMULTIPLEXING_BARCODE_LOADER_HH

#include <boost/mpl/equal_to.hpp>

#include "common/Debug.hh"
#include "common/Threads.hpp"

#include "demultiplexing/Barcode.hh"

#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/TileMetadata.hh"

#include "io/BclMapper.hh"
#include "io/FileBufCache.hh"

namespace isaac
{
namespace demultiplexing
{

/**
 ** \brief Encapsulates the variables that are shared by all the threads while
 ** loading the barcodes.
 **/
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
        const bool ignoreMissingBcls,
        const flowcell::TileMetadataList &tileMetadataList,
        const std::vector<unsigned> &barcodeCycles,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const unsigned maxThreads);

    void load(std::vector<Barcode> &barcodes,
              std::vector<Barcode>::iterator &nextTileBarcodes,
              flowcell::TileMetadataList::const_iterator &nextTile,
              const flowcell::TileMetadataList::const_iterator tilesEnd,
              const unsigned threadNumber);
private:
    // The mutex used to acquire the next tile and the destination of the seeds
    boost::mutex mutex_;
    const flowcell::TileMetadataList &tileMetadataList_;
    const std::vector<unsigned> barcodeCycles_;
    const unsigned unknownBarcodeIndex_;

    boost::ptr_vector<io::SingleCycleBclMapper> threadBclMappers_;

    void loadTileCycle(
        io::SingleCycleBclMapper &bclMapper,
        std::vector<Barcode>::iterator destination,
        const flowcell::TileMetadata &tile,
        const unsigned cycle);
};

class BarcodeLoader: boost::noncopyable
{
public:
    BarcodeLoader(
        const bool ignoreMissingBcls,
        common::ThreadVector &threads,
        const unsigned inputLoadersMax,
        const flowcell::TileMetadataList &allTilesMetadata,
        const std::vector<unsigned> &barcodeCycles,
        const flowcell::BarcodeMetadataList &barcodeMetadataList);

    static void allocate(const flowcell::TileMetadataList &tiles, std::vector<Barcode> &barcodes);
    static bool seeIfFits(const flowcell::TileMetadataList &tiles);
    static bool selectTiles(
        flowcell::TileMetadataList &unprocessedPool,
        flowcell::TileMetadataList &selectedTiles);

    void loadBarcodes(const flowcell::TileMetadataList &tiles, std::vector<Barcode> &barcodes);

private:
    const unsigned inputLoadersMax_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;

    common::ThreadVector &threads_;
    ParallelBarcodeLoader parallelBarcodeLoader_;
};

} // namespace demultiplexing
} // namespace isaac

#endif // #ifndef iSAAC_DEMULTIPLEXING_BARCODE_LOADER_HH
