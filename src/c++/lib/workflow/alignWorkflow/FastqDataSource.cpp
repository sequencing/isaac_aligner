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
 ** \file FastqDataSource.cpp
 **
 ** \brief see FastqDataSource.hh
 **
 ** \author Roman Petrovski
 **/

#include "FastqDataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

unsigned FastqDataSource::determineMemoryCapacity(const unsigned clusterLength)
{
    ISAAC_THREAD_CERR << "Determining memory capacity for Fastq data" << std::endl;
    // Allocate as much memory as possible in tileClustersMax_ decrements
    alignment::BclClusters testClusters(clusterLength);
    // Assume seeds will equate at most to the cluster length at a time. This
    // means the RAM available for clusters will be third (seeds are made of kmer and metadata)
    // the ram successfully allocated
    // TODO: start with something more reasonable than hard-cored 256 Gigs
    size_t testCount = 1024UL * 1024UL * 1024UL * 256UL / clusterLength;
    while (testCount)
    {
        try
        {
            testClusters.reserveClusters(testCount);
            break;
        }
        catch (std::bad_alloc &e)
        {
            if (tileClustersMax_ >= testCount)
            {
                BOOST_THROW_EXCEPTION(common::MemoryException((boost::format(
                    "Insufficient memory to load seeds even for %u fastq clusters") % testCount).str()));
            }
            errno = 0;
            testCount -= tileClustersMax_;
        }
    }

    ISAAC_THREAD_CERR << "Determining memory capacity for Fastq data done. " <<
        testCount << " clusters of length " << clusterLength <<
            " will fit" << std::endl;
    return testCount;
}

FastqDataSource::FastqDataSource(
    const bool allowVariableLength,
    const unsigned coresMax,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const reference::SortedReferenceXmlList &sortedReferenceXmlList,
    const flowcell::Layout &fastqFlowcellLayout) :
        coresMax_(coresMax),
        barcodeMetadataList_(barcodeMetadataList),
        fastqFlowcellLayout_(fastqFlowcellLayout),
        sortedReferenceXmlList_(sortedReferenceXmlList),
        clusterLength_(flowcell::getTotalReadLength(fastqFlowcellLayout.getReadMetadataList())),
        clustersAtATimeMax_(determineMemoryCapacity(clusterLength_)),
        clusters_(clusterLength_),
        lanes_(fastqFlowcellLayout.getLaneIds()),
        currentLaneIterator_(lanes_.begin()),
        currentTile_(1),
        fastqLoader_(allowVariableLength)

{
}


flowcell::TileMetadataList FastqDataSource::discoverTiles()
{
    loadedTiles_.clear();

    const bool compressed = fastqFlowcellLayout_.getFormat() == flowcell::Layout::FastqGz;

    // As we don't know whether the current lane has been completely loaded or we're in the
    // middle of discovering it's tiles, just attempt to load more data for it and stop only
    // when some data is loaded for this or subsequent lane or we run out of lanes.
    unsigned clustersLoaded = 0;
    while (!clustersLoaded)
    {
        if (lanes_.end() == currentLaneIterator_)
        {
            return loadedTiles_;
        }

        // Allocate one third the available RAM for clusters as for each cluster kmer we need a kmer+metadata for the seed
        const unsigned clustersToLoad = clustersAtATimeMax_ / 3;
        clusters_.reset(clusterLength_, clustersToLoad);
        // load clusters, return tile breakdown based on tileClustersMax_

        boost::filesystem::path read1Path;
        flowcell::Layout::getFastqFilePath(1, *currentLaneIterator_, fastqFlowcellLayout_.getBaseCallsDirectory(),
                                           compressed, read1Path);
        boost::filesystem::path read2Path;
        flowcell::Layout::getFastqFilePath(2, *currentLaneIterator_, fastqFlowcellLayout_.getBaseCallsDirectory(),
                                           compressed, read2Path);
        // this will keep the current files open if the paths don't change
        fastqLoader_.open(read1Path, read2Path);
        std::vector<char>::iterator clustersEnd = clusters_.cluster(0);
        clustersLoaded = fastqLoader_.loadClusters(clustersToLoad, fastqFlowcellLayout_.getReadMetadataList(), clustersEnd);
        ISAAC_THREAD_CERR << "Loaded  " << clustersLoaded << " clusters of length " << clusterLength_ << std::endl;
        clusters_.reset(clusterLength_, clustersLoaded);
        ISAAC_ASSERT_MSG(clusters_.end() == clustersEnd, "Mismatch between expected and actual cluster buffer ends. Distance: "<<
                         std::distance(clustersEnd, clusters_.end()));

        if (!clustersLoaded)
        {
            // there was nothing left to load for this lane, proceed with the next one.
            ++currentLaneIterator_;
            currentTile_ = 1;
        }
    }

    const flowcell::TileMetadata::Compression compression =
        compressed ? flowcell::TileMetadata::GzCompression : flowcell::TileMetadata::NoCompression;

    const std::string &flowcellId = fastqFlowcellLayout_.getFlowcellId();
    while (true)
    {
        const unsigned tileClustersMax = tileClustersMax_; //otherwise linker fails with gcc 4.6.1 Debug builds
        const unsigned clusterCount = std::min(clustersLoaded, tileClustersMax);
        const flowcell::TileMetadata tileMetadata(
            flowcellId, fastqFlowcellLayout_.getIndex(),
            currentTile_++, *currentLaneIterator_, fastqFlowcellLayout_.getBaseCallsDirectory(),
            clusterCount,
            compression,
            loadedTiles_.size());
        loadedTiles_.push_back(tileMetadata);

        if (clustersLoaded < tileClustersMax_)
        {
            break;
        }
        clustersLoaded -= tileClustersMax_;
    }

    return loadedTiles_;
}

void FastqDataSource::initBuffers(
    flowcell::TileMetadataList &unprocessedTiles,
    const alignment::SeedMetadataList &seedMetadataList,
    common::ThreadVector &threads)
{
    seedGenerator_.reset(new alignment::ClusterSeedGenerator(
        threads,
        coresMax_, barcodeMetadataList_,
        fastqFlowcellLayout_.getReadMetadataList(),
        seedMetadataList,
        sortedReferenceXmlList_, unprocessedTiles,
        clusters_,
        loadedTiles_));
}

void FastqDataSource::generateSeeds(
    const flowcell::TileMetadataList &tiles,
    const alignment::matchFinder::TileClusterInfo &tileClusterBarcode,
    std::vector<alignment::Seed> &seeds,
    common::ScoopedMallocBlock  &mallocBlock)
{
    seedGenerator_->generateSeeds(tiles, tileClusterBarcode, seeds, mallocBlock);
}

const std::vector<std::vector<alignment::Seed>::iterator> &FastqDataSource::getReferenceSeedBounds()
{
    return seedGenerator_->getReferenceSeedBounds();
}

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
