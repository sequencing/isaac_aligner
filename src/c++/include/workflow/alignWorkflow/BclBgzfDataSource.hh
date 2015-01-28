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
 ** \file BclBgzfDataSource.hh
 **
 ** \brief Encapsulation of BaseCalls folder with bcl.bgzf and filter files as seed and cluster data source
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCLBGZF_DATA_SOURCE_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCLBGZF_DATA_SOURCE_HH

#include "alignment/matchFinder/TileClusterInfo.hh"
#include "alignment/SeedLoader.hh"
#include "demultiplexing/BarcodeLoader.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/TileMetadata.hh"
#include "io/LocsMapper.hh"
#include "io/ClocsMapper.hh"
#include "io/FiltersMapper.hh"
#include "rta/BclBgzfTileReader.hh"
#include "rta/BclMapper.hh"
#include "rta/LaneBciMapper.hh"
#include "workflow/alignWorkflow/DataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

class BclBgzfTileSource
{
public:

    unsigned readTileCycle(
        const flowcell::TileMetadata &tile,
        char *cycleBuffer,
        const std::size_t bufferSize) const;
};

template <typename KmerT>
class BclBgzfSeedSource : public TileSource, public BarcodeSource, public SeedSource<KmerT>
{
    typedef alignment::Seed<KmerT> SeedT;
    typedef typename std::vector<SeedT>::iterator SeedIterator;

    const bool ignoreMissingBcls_;
    const unsigned inputLoadersMax_;
    const unsigned coresMax_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::Layout &bclFlowcellLayout_;
    common::ThreadVector &threads_;
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList_;
    /// mapping from tile index to the index of the tile in bci.
    std::vector<unsigned> tileBciIndexMap_;
    const flowcell::TileMetadataList flowcellTiles_;
    const unsigned maxTileClusterCount_;
    boost::scoped_ptr<alignment::ParallelSeedLoader<rta::BclBgzfTileReader, KmerT> > seedLoader_;
    // holds the state across multiple discoverTiles calls
    flowcell::TileMetadataList::const_iterator undiscoveredTiles_;
    std::vector<rta::BclBgzfTileReader> threadBclReaders_;
    demultiplexing::BarcodeLoader<rta::BclBgzfTileReader> barcodeLoader_;

    boost::ptr_vector<rta::SingleCycleBclMapper<rta::BclBgzfTileReader> > threadBclMappers_;
    // geometry [cycle]
    std::vector<rta::CycleBciMapper> cycleBciMappers_;
public:
    BclBgzfSeedSource(
        const bool ignoreMissingBcls,
        const unsigned inputLoadersMax,
        const unsigned coresMax,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const flowcell::Layout &bclFlowcellLayout,
        common::ThreadVector &threads);

    // TileSource implementation
    flowcell::TileMetadataList discoverTiles();

    // BarcodeSource implementation
    virtual void loadBarcodes(
        const unsigned unknownBarcodeIndex,
        const flowcell::TileMetadataList &tiles,
        std::vector<demultiplexing::Barcode> &barcodes);

    // SeedSource implementation
    void initBuffers(
        flowcell::TileMetadataList &unprocessedTiles,
        const alignment::SeedMetadataList &seedMetadataList);
    void generateSeeds(
        const flowcell::TileMetadataList &tiles,
        const alignment::matchFinder::TileClusterInfo &tileClusterBarcode,
        std::vector<SeedT > &seeds,
        common::ScoopedMallocBlock  &mallocBlock);
    const std::vector<SeedIterator> &getReferenceSeedBounds() const;
private:
    /**
     * \return vector of tiles ordered by: flowcellId_, lane_, tile_
     */
    static flowcell::TileMetadataList getTiles(
        const flowcell::Layout &flowcellLayout,
        std::vector<unsigned> &tileBciIndexMap);
    void initCycleBciMappers(
        const std::vector<unsigned>& cycles,
        const unsigned laneNumber,
        std::vector<rta::CycleBciMapper> &cycleBciMappers) const;
};


class BclBgzfBaseCallsSource : boost::noncopyable
{
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    common::ThreadVector &bclLoadThreads_;

    // temporaries to avoid memory allocations during data processing
    boost::filesystem::path bciFilePath_;
    boost::filesystem::path filterFilePath_;
    boost::filesystem::path positionsFilePath_;
    std::vector<unsigned> cycles_;
    io::FileBufWithReopen cycleBciFileBuf_;

    // References to cycleBciMappers_ and tileBciIndexMap_ get passed to threadReaders_ during construction.
    // The contents is dynamically swapped each time we switch flowcell lane
    /// Cumulative offsets of each tile in clusters. All tiles of the lane included
    std::vector<unsigned long> tileClusterOffsets_;
    std::vector<rta::CycleBciMapper> cycleBciMappers_;
    std::vector<unsigned> tileBciIndexMap_;

    rta::LaneBciMapper laneBciMapper_;

    std::vector<rta::BclBgzfTileReader> threadReaders_;
    rta::ParallelBclMapper<rta::BclBgzfTileReader> bclMapper_;
    io::FiltersMapper filtersMapper_;
    io::ClocsMapper clocsMapper_;
    io::LocsMapper locsMapper_;


    unsigned currentFlowcellIndex_;
    unsigned currentLaneNumber_;



public:
    BclBgzfBaseCallsSource(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const flowcell::TileMetadataList &tileMetadataList,
        const bool ignoreMissingBcls,
        const bool ignoreMissingFilters,
        common::ThreadVector &bclLoadThreads,
        const unsigned inputLoadersMax,
        const bool extractClusterXy);
    void loadClusters(
        const flowcell::TileMetadataList &allTiles,
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData);

private:
    void initBciMappers(
        const flowcell::Layout &flowcell,
        const flowcell::TileMetadataList &allTiles,
        const unsigned laneNumber,
        std::vector<rta::CycleBciMapper>  &cycleBciMappers,
        std::vector<unsigned long> &tileClusterOffsets,
        std::vector<unsigned> &tileBciIndexMap);

    void bclToClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData,
        const bool useLocsPositions) const;
};


} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCLBGZF_DATA_SOURCE_HH
