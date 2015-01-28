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
 ** \file FastqDataSource.hh
 **
 ** \brief Encapsulation of single-ended and paired data stored in fastq file(s)
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FASTQ_DATA_SOURCE_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_FASTQ_DATA_SOURCE_HH

#include "alignment/matchFinder/TileClusterInfo.hh"
#include "alignment/BclClusters.hh"
#include "alignment/ClusterSeedGenerator.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/FastqLayout.hh"
#include "flowcell/TileMetadata.hh"
#include "io/FastqLoader.hh"
#include "workflow/alignWorkflow/DataSource.hh"


namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

template <typename KmerT>
class FastqSeedSource : public TileSource, public BarcodeSource, public SeedSource<KmerT>
{
    typedef alignment::Seed<KmerT> SeedT;
    typedef typename std::vector<SeedT>::iterator SeedIterator;

    const unsigned tileClustersMax_;
    const unsigned coresMax_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::Layout &fastqFlowcellLayout_;
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList_;
    const unsigned clusterLength_;
    const unsigned clustersAtATimeMax_;
    alignment::BclClusters clusters_;
    flowcell::TileMetadataList loadedTiles_;
    const std::vector<unsigned> lanes_;
    std::vector<unsigned>::const_iterator currentLaneIterator_;
    unsigned currentTile_;
    common::ThreadVector &threads_;
    io::FastqLoader fastqLoader_;
    boost::scoped_ptr<alignment::ClusterSeedGenerator<KmerT> > seedGenerator_;

public:
    FastqSeedSource(
        const unsigned long availableMemory,
        const unsigned clustersAtATimeMax,
        const bool allowVariableLength,
        const unsigned coresMax,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const flowcell::Layout &fastqFlowcellLayout,
        common::ThreadVector &threads);

    // TileSource implementation
    flowcell::TileMetadataList discoverTiles();

    // BarcodeSource implementation
    virtual void loadBarcodes(
        const unsigned unknownBarcodeIndex,
        const flowcell::TileMetadataList &tiles,
        std::vector<demultiplexing::Barcode> &barcodes)
    {
        ISAAC_ASSERT_MSG(false, "Barcode resolution is not implemented for Fastq data");
    }

    // SeedSource implementation
    void initBuffers(
        flowcell::TileMetadataList &unprocessedTiles,
        const alignment::SeedMetadataList &seedMetadataList);
    void generateSeeds(
        const flowcell::TileMetadataList &tiles,
        const alignment::matchFinder::TileClusterInfo &tileClusterBarcode,
        std::vector<SeedT> &seeds,
        common::ScoopedMallocBlock  &mallocBlock);
    const std::vector<SeedIterator> &getReferenceSeedBounds() const;

private:
    static unsigned determineMemoryCapacity(
        const unsigned long availableMemory,
        const unsigned tileClustersMax,
        const std::size_t seedsPerCluster,
        const unsigned clusterLength);
};


class FastqBaseCallsSource : boost::noncopyable
{
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    io::FastqLoader fastqLoader_;
    std::vector<boost::filesystem::path> fastqFilePaths_;

public:
    FastqBaseCallsSource(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const flowcell::TileMetadataList &tileMetadataList,
        const bool allowVariableFastqLength,
        common::ThreadVector &threads,
        const unsigned inputLoadersMax);
    void loadClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData);
};

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FASTQ_DATA_SOURCE_HH
