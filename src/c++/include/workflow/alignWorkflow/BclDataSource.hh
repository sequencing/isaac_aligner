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
 ** \file BclDataSource.hh
 **
 ** \brief Encapsulation of BaseCalls folder with bcl and filter files as seed and cluster data source
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCL_DATA_SOURCE_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCL_DATA_SOURCE_HH

#include "alignment/matchFinder/TileClusterInfo.hh"
#include "alignment/SeedLoader.hh"
#include "demultiplexing/BarcodeLoader.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/TileMetadata.hh"
#include "io/LocsMapper.hh"
#include "io/ClocsMapper.hh"
#include "io/FiltersMapper.hh"
#include "rta/BclReader.hh"
#include "rta/BclMapper.hh"
#include "workflow/alignWorkflow/DataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

template <typename KmerT>
class BclSeedSource : public TileSource, public BarcodeSource, public SeedSource<KmerT>
{
    typedef alignment::Seed<KmerT> SeedT;
    typedef typename std::vector<SeedT>::iterator SeedIterator;

    const bool ignoreMissingBcls_;
    const unsigned inputLoadersMax_;
    const unsigned coresMax_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::Layout &bclFlowcellLayout_;
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList_;
    const flowcell::TileMetadataList flowcellTiles_;
    const unsigned maxTileClusterCount_;
    boost::scoped_ptr<alignment::ParallelSeedLoader<rta::BclReader, KmerT> > seedLoader_;
    // holds the state across multiple discoverTiles calls
    flowcell::TileMetadataList::const_iterator undiscoveredTiles_;
    std::vector<rta::BclReader> threadBclReaders_;
    common::ThreadVector &threads_;
    const unsigned longestBclPathLength_;
    demultiplexing::BarcodeLoader<rta::BclReader> barcodeLoader_;

    boost::ptr_vector<rta::SingleCycleBclMapper<rta::BclReader> > threadBclMappers_;
public:
    BclSeedSource(
        const bool ignoreMissingBcls,
        const unsigned inputLoadersMax,
        const unsigned coresMax,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const flowcell::Layout &bclFlowcellLayout,
        common::ThreadVector &threads);

    // TileSource implementation
    virtual flowcell::TileMetadataList discoverTiles();

    // BarcodeSource implementation
    virtual void loadBarcodes(
        const unsigned unknownBarcodeIndex,
        const flowcell::TileMetadataList &tiles,
        std::vector<demultiplexing::Barcode> &barcodes);

    // SeedSource implementation
    virtual void initBuffers(
        flowcell::TileMetadataList &unprocessedTiles,
        const alignment::SeedMetadataList &seedMetadataList);
    virtual  void generateSeeds(
        const flowcell::TileMetadataList &tiles,
        const alignment::matchFinder::TileClusterInfo &tileClusterBarcode,
        std::vector<SeedT > &seeds,
        common::ScoopedMallocBlock  &mallocBlock);
    const std::vector<SeedIterator> &getReferenceSeedBounds() const;
private:
    /**
     * \return vector of tiles ordered by: flowcellId_, lane_, tile_
     */
    flowcell::TileMetadataList getTiles(const flowcell::Layout &flowcellLayout) const;

};


class BclBaseCallsSource : boost::noncopyable
{
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    common::ThreadVector &bclLoadThreads_;
    boost::filesystem::path filterFilePath_;
    boost::filesystem::path positionsFilePath_;
    std::vector<rta::BclReader> threadReaders_;
    rta::ParallelBclMapper<rta::BclReader> bclMapper_;
    io::FiltersMapper filtersMapper_;
    io::ClocsMapper clocsMapper_;
    io::LocsMapper locsMapper_;

public:
    BclBaseCallsSource(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const flowcell::TileMetadataList &tileMetadataList,
        const bool ignoreMissingBcls,
        const bool ignoreMissingFilters,
        common::ThreadVector &bclLoadThreads,
        const unsigned inputLoadersMax,
        const bool extractClusterXy);
    void loadClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData);

private:
    void bclToClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData,
        const bool useLocsPositions) const;
};


} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCL_DATA_SOURCE_HH
