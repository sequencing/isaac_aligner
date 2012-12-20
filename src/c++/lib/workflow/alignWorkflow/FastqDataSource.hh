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
#include "flowcell/Layout.hh"
#include "flowcell/TileMetadata.hh"
#include "io/FastqLoader.hh"
#include "workflow/alignWorkflow/DataSource.hh"


namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

class FastqDataSource : public DataSource
{
    static const unsigned tileClustersMax_ = 4000000;
    const unsigned coresMax_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::Layout &fastqFlowcellLayout_;
    const reference::SortedReferenceXmlList &sortedReferenceXmlList_;
    const unsigned clusterLength_;
    const unsigned clustersAtATimeMax_;
    alignment::BclClusters clusters_;
    flowcell::TileMetadataList loadedTiles_;
    const std::vector<unsigned> lanes_;
    std::vector<unsigned>::const_iterator currentLaneIterator_;
    unsigned currentTile_;
    io::FastqLoader fastqLoader_;
    boost::scoped_ptr<alignment::ClusterSeedGenerator> seedGenerator_;

public:
    FastqDataSource(
        const bool allowVariableLength,
        const unsigned coresMax,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const reference::SortedReferenceXmlList &sortedReferenceXmlList,
        const flowcell::Layout &fastqFlowcellLayout);

    // DataSource implementation

    flowcell::TileMetadataList discoverTiles();
    void initBuffers(
        flowcell::TileMetadataList &unprocessedTiles,
        const alignment::SeedMetadataList &seedMetadataList,
        common::ThreadVector &threads);
    void generateSeeds(
        const flowcell::TileMetadataList &tiles,
        const alignment::matchFinder::TileClusterInfo &tileClusterBarcode,
        std::vector<alignment::Seed> &seeds,
        common::ScoopedMallocBlock  &mallocBlock);
    const std::vector<std::vector<alignment::Seed>::iterator> &getReferenceSeedBounds();

private:
    static unsigned determineMemoryCapacity(const unsigned clusterLength);
};

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FASTQ_DATA_SOURCE_HH
