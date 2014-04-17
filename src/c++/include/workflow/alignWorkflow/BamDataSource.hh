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
 ** \file BamDataSource.hh
 **
 ** \brief Encapsulation of single-ended and paired data stored in bam file
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BAM_DATA_SOURCE_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_BAM_DATA_SOURCE_HH

#include "alignment/matchFinder/TileClusterInfo.hh"
#include "alignment/BclClusters.hh"
#include "alignment/ClusterSeedGenerator.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/BamLayout.hh"
#include "flowcell/TileMetadata.hh"
#include "io/BamLoader.hh"
#include "workflow/alignWorkflow/bamDataSource/PairedEndClusterExtractor.hh"
#include "workflow/alignWorkflow/DataSource.hh"


namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

class BamClusterLoader
{
    std::string flowcellId_;

    io::BamLoader bamLoader_;
    bamDataSource::PairedEndClusterExtractor clusterExtractor_;
public:
    BamClusterLoader(
        const bool cleanupIntermediary,
        const std::size_t maxPathLength,
        common::ThreadVector &threads,
        const unsigned coresMax,
        const boost::filesystem::path &tempDirectoryPath,
        const std::size_t maxBamFileLength,
        const std::size_t maxFlowcellIdLength,
        const std::size_t minClusterLength) :
        bamLoader_(maxPathLength, threads, coresMax),
        clusterExtractor_(tempDirectoryPath, maxBamFileLength, maxFlowcellIdLength, minClusterLength, cleanupIntermediary)
    {
        flowcellId_.reserve(maxFlowcellIdLength);
    }

    void open(
        const std::string &flowcellId,
        const boost::filesystem::path &bamPath);

    template <typename ClusterInsertIt, typename PfInserIt>
    unsigned loadClusters(unsigned clusterCount, const flowcell::ReadMetadataList &readMetadataList,
        ClusterInsertIt &clusterIt, PfInserIt &pfIt);

private:
    template <typename InsertIt, typename PfInserIt>
    unsigned loadPairedReads(unsigned clusterCount, const flowcell::ReadMetadataList &readMetadataList, InsertIt &it, PfInserIt &pfIt);
    template <typename InsertIt, typename PfInserIt>
    unsigned loadSingleReads(unsigned clusterCount, const flowcell::ReadMetadataList &readMetadataList, InsertIt &it, PfInserIt &pfIt);
};


template <typename KmerT>
class BamSeedSource : public TileSource, public BarcodeSource, public SeedSource<KmerT>
{
    typedef alignment::Seed<KmerT> SeedT;
    typedef typename std::vector<SeedT>::iterator SeedIterator;

    const flowcell::Layout &bamFlowcellLayout_;
    const unsigned tileClustersMax_;
    const unsigned coresMax_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList_;
    const unsigned clusterLength_;
    const unsigned clustersAtATimeMax_;
    alignment::BclClusters clusters_;
    flowcell::TileMetadataList loadedTiles_;
    unsigned currentTile_;
    common::ThreadVector &threads_;
    BamClusterLoader bamClusterLoader_;
    boost::scoped_ptr<alignment::ClusterSeedGenerator<KmerT> > seedGenerator_;

public:
    BamSeedSource(
        const boost::filesystem::path &tempDirectoryPath,
        const unsigned long availableMemory,
        const unsigned clustersAtATimeMax,
        const bool cleanupIntermediary,
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
        ISAAC_ASSERT_MSG(false, "Barcode resolution is not implemented for Bam data");
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
        const unsigned clusterLength);
};


class BamBaseCallsSource : boost::noncopyable
{
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    BamClusterLoader bamClusterLoader_;
    boost::filesystem::path bamFilePath_;

public:
    BamBaseCallsSource(
        const boost::filesystem::path &tempDirectoryPath,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const flowcell::TileMetadataList &tileMetadataList,
        const bool cleanupIntermediary,
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
