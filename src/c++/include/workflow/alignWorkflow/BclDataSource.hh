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
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/TileMetadata.hh"
#include "io/BclMapper.hh"
#include "io/ClocsMapper.hh"
#include "io/FiltersMapper.hh"
#include "workflow/alignWorkflow/DataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

template <typename KmerT>
class BclSeedSource : public SeedSource<KmerT>
{
    typedef alignment::Seed<KmerT> SeedT;
    typedef typename std::vector<SeedT>::iterator SeedIterator;

    const bool ignoreMissingBcls_;
    const unsigned inputLoadersMax_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::Layout &bclFlowcellLayout_;
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList_;
    const flowcell::TileMetadataList flowcellTiles_;
    boost::scoped_ptr<alignment::ParallelSeedLoader<KmerT> > seedLoader_;
    // holds the state across multiple discoverTiles calls
    flowcell::TileMetadataList::const_iterator undiscoveredTiles_;
public:
    BclSeedSource(
        const bool ignoreMissingBcls,
        const unsigned inputLoadersMax,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const flowcell::Layout &bclFlowcellLayout) :
            ignoreMissingBcls_(ignoreMissingBcls),
            inputLoadersMax_(inputLoadersMax),
            barcodeMetadataList_(barcodeMetadataList),
            bclFlowcellLayout_(bclFlowcellLayout),
            sortedReferenceMetadataList_(sortedReferenceMetadataList),
            flowcellTiles_(getTiles(bclFlowcellLayout)),
            undiscoveredTiles_(flowcellTiles_.begin())
    {
    }

    // SeedSource implementation

    flowcell::TileMetadataList discoverTiles();
    void initBuffers(
        flowcell::TileMetadataList &unprocessedTiles,
        const alignment::SeedMetadataList &seedMetadataList,
        common::ThreadVector &threads);
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
    flowcell::TileMetadataList getTiles(const flowcell::Layout &flowcellLayout) const;

};


class BclBaseCallsSource : boost::noncopyable
{
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    common::ThreadVector &bclLoadThreads_;
    std::vector<boost::filesystem::path> bclFilePaths_;
    boost::filesystem::path filterFilePath_;
    boost::filesystem::path positionsFilePath_;
    io::ParallelBclMapper bclMapper_;
    io::FiltersMapper filtersMapper_;
    io::ClocsMapper positionsMapper_;

public:
    BclBaseCallsSource(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const flowcell::TileMetadataList &tileMetadataList,
        const bool ignoreMissingBcls,
        const bool ignoreMissingFilters,
        common::ThreadVector &bclLoadThreads,
        const unsigned inputLoadersMax);
    void loadClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData);

private:
    void bclToClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData) const;
};


} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCL_DATA_SOURCE_HH
