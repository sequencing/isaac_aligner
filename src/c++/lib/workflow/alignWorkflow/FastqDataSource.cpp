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
 ** \file FastqDataSource.cpp
 **
 ** \brief see FastqDataSource.hh
 **
 ** \author Roman Petrovski
 **/

#include "workflow/alignWorkflow/FastqDataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

template <typename KmerT>
unsigned FastqSeedSource<KmerT>::determineMemoryCapacity(
    const unsigned long availableMemory,
    const unsigned tileClustersMax,
    const std::size_t seedsPerCluster,
    const unsigned clusterLength)
{
    ISAAC_THREAD_CERR << "Determining memory capacity for Fastq data" << std::endl;
    // Allocate as much memory as possible in tileClustersMax_ decrements
    alignment::BclClusters testClusters(clusterLength);
    // Assume seeds will equate at most to the cluster length at a time. This
    // means the RAM available for clusters will be one third (seeds are made of kmer and metadata)
    // the ram successfully allocated
    size_t testCount = availableMemory / clusterLength;
    while (testCount)
    {
        try
        {
            std::vector<alignment::Seed<KmerT> > seeds;
            seeds.reserve(testCount * seedsPerCluster);
//            ISAAC_THREAD_CERR << "Determining memory capacity trying. " << testCount << std::endl;
            testClusters.reserveClusters(testCount, false);
            break;
        }
        catch (std::bad_alloc &e)
        {
            if (tileClustersMax >= testCount)
            {
                BOOST_THROW_EXCEPTION(common::MemoryException((boost::format(
                    "Insufficient memory to load seeds even for %u fastq clusters") % testCount).str()));
            }
            // reset errno, to prevent misleading error messages when failing code does not set errno
            errno = 0;
            testCount -= tileClustersMax;
        }
    }

    ISAAC_THREAD_CERR << "Determining memory capacity for Fastq data done. " <<
        testCount << " clusters of length " << clusterLength <<
            " will fit" << std::endl;
    return testCount;
}

template <typename KmerT>
FastqSeedSource<KmerT>::FastqSeedSource(
    const unsigned long availableMemory,
    const unsigned clustersAtATimeMax,
    const bool allowVariableLength,
    const unsigned coresMax,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
    const flowcell::Layout &fastqFlowcellLayout,
    common::ThreadVector &threads) :
        tileClustersMax_(clustersAtATimeMax ?
            std::min<unsigned>(clustersAtATimeMax, 40000000 / fastqFlowcellLayout.getSeedMetadataList().size()):
            (40000000 / fastqFlowcellLayout.getSeedMetadataList().size())),
        coresMax_(coresMax),
        barcodeMetadataList_(barcodeMetadataList),
        fastqFlowcellLayout_(fastqFlowcellLayout),
        sortedReferenceMetadataList_(sortedReferenceMetadataList),
        clusterLength_(flowcell::getTotalReadLength(fastqFlowcellLayout.getReadMetadataList())),
        clustersAtATimeMax_(clustersAtATimeMax ? clustersAtATimeMax : determineMemoryCapacity(availableMemory, tileClustersMax_, fastqFlowcellLayout.getSeedMetadataList().size(), clusterLength_)),
        clusters_(clusterLength_),
        lanes_(fastqFlowcellLayout.getLaneIds()),
        currentLaneIterator_(lanes_.begin()),
        currentTile_(1),
        threads_(threads),
        fastqLoader_(allowVariableLength, 0, threads_, coresMax_)

{
}


template <typename KmerT>
flowcell::TileMetadataList FastqSeedSource<KmerT>::discoverTiles()
{
    loadedTiles_.clear();

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

        const unsigned clustersToLoad = clustersAtATimeMax_;
        clusters_.reset(clusterLength_, clustersToLoad);
        // load clusters, return tile breakdown based on tileClustersMax_

        const boost::filesystem::path read1Path =
            fastqFlowcellLayout_.getLaneReadAttribute<flowcell::Layout::Fastq, flowcell::FastqFilePathAttributeTag>(
                *currentLaneIterator_, fastqFlowcellLayout_.getReadMetadataList().at(0).getNumber());
        if (1 == fastqFlowcellLayout_.getReadMetadataList().size())
        {
            // this will keep the current files open if the paths don't change
            fastqLoader_.open(read1Path);
        }
        else // assume paired data
        {
            const boost::filesystem::path read2Path =
                fastqFlowcellLayout_.getLaneReadAttribute<flowcell::Layout::Fastq, flowcell::FastqFilePathAttributeTag>(
                    *currentLaneIterator_, fastqFlowcellLayout_.getReadMetadataList().at(1).getNumber());
            // this will keep the current files open if the paths don't change
            fastqLoader_.open(read1Path, read2Path);
        }

        std::vector<char>::iterator clustersEnd = clusters_.cluster(0);
        clustersLoaded = fastqLoader_.loadClusters(clustersToLoad, fastqFlowcellLayout_.getReadMetadataList(), clustersEnd);
        ISAAC_THREAD_CERR << "Loaded  " << clustersLoaded << " clusters of length " << clusterLength_ << std::endl;
        clusters_.reset(clusterLength_, clustersLoaded);
        if (clustersLoaded < clustersToLoad / 2)
        {
            //allocated too much memory for bcl data. Might pay off to give some away especially if seeds are placed very densely for high sensitivity
            clusters_.reduceWastedMemory();
        }

        if (!clustersLoaded)
        {
            // there was nothing left to load for this lane, proceed with the next one.
            ++currentLaneIterator_;
            currentTile_ = 1;
        }
    }

    const std::string &flowcellId = fastqFlowcellLayout_.getFlowcellId();
    while (clustersLoaded)
    {
        const unsigned clusterCount = std::min(clustersLoaded, tileClustersMax_);
        const flowcell::TileMetadata tileMetadata(
            flowcellId, fastqFlowcellLayout_.getIndex(),
            currentTile_++, *currentLaneIterator_,
            clusterCount,
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

template <typename KmerT>
void FastqSeedSource<KmerT>::initBuffers(
    flowcell::TileMetadataList &unprocessedTiles,
    const alignment::SeedMetadataList &seedMetadataList)
{
    seedGenerator_.reset(new alignment::ClusterSeedGenerator<KmerT>(
        threads_,
        coresMax_, barcodeMetadataList_,
        fastqFlowcellLayout_,
        seedMetadataList,
        sortedReferenceMetadataList_, unprocessedTiles,
        clusters_,
        loadedTiles_));
}

template <typename KmerT>
void FastqSeedSource<KmerT>::generateSeeds(
    const flowcell::TileMetadataList &tiles,
    const alignment::matchFinder::TileClusterInfo &tileClusterBarcode,
    std::vector<SeedT> &seeds,
    common::ScoopedMallocBlock  &mallocBlock)
{
    seedGenerator_->generateSeeds(tiles, tileClusterBarcode, seeds, mallocBlock);
}

template <typename KmerT>
const std::vector<typename FastqSeedSource<KmerT>::SeedIterator> &FastqSeedSource<KmerT>::getReferenceSeedBounds() const
{
    return seedGenerator_->getReferenceSeedBounds();
}

/////////////// FastqBaseCallsSource implementation
inline boost::filesystem::path getLongestFastqPath(const flowcell::FlowcellLayoutList &flowcellLayoutList)
{
    return flowcell::getLongestAttribute<flowcell::Layout::Fastq, flowcell::FastqFilePathAttributeTag>(flowcellLayoutList);
}

FastqBaseCallsSource::FastqBaseCallsSource(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const flowcell::TileMetadataList &tileMetadataList,
    const bool allowVariableFastqLength,
    common::ThreadVector &threads,
    const unsigned inputLoadersMax):
    flowcellLayoutList_(flowcellLayoutList),
    fastqLoader_(allowVariableFastqLength, getLongestFastqPath(flowcellLayoutList_).string().size(), threads, inputLoadersMax),
    fastqFilePaths_(2) //read 1 and read 2 paths
{
    boost::filesystem::path longestFastqFilePath = getLongestFastqPath(flowcellLayoutList_);
    BOOST_FOREACH(boost::filesystem::path &p, fastqFilePaths_)
    {
        // this has to be done separately for each path or else they all share one buffer
        p = longestFastqFilePath.c_str();
    }
}

void FastqBaseCallsSource::loadClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData)
{
    ISAAC_THREAD_CERR << "Loading Fastq data for " << tileMetadata << std::endl;
    const flowcell::Layout &flowcell = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex());

    flowcell.getLaneReadAttribute<flowcell::Layout::Fastq, flowcell::FastqFilePathAttributeTag>(
        tileMetadata.getLane(), flowcell.getReadMetadataList().at(0).getNumber(), fastqFilePaths_.at(0));

    if (1 == flowcell.getReadMetadataList().size())
    {
        fastqLoader_.open(fastqFilePaths_[0]);
    }
    else
    {
        flowcell.getLaneReadAttribute<flowcell::Layout::Fastq, flowcell::FastqFilePathAttributeTag>(
            tileMetadata.getLane(), flowcell.getReadMetadataList().at(1).getNumber(), fastqFilePaths_.at(1));
        fastqLoader_.open(fastqFilePaths_[0], fastqFilePaths_[1]);
    }

    const unsigned clustersToLoad = tileMetadata.getClusterCount();
    ISAAC_THREAD_CERR << "Resetting Fastq data for " << clustersToLoad << " clusters" << std::endl;
    bclData.reset(flowcell::getTotalReadLength(flowcell.getReadMetadataList()), clustersToLoad);
    ISAAC_THREAD_CERR << "Resetting Fastq data done for " << bclData.getClusterCount() << " clusters" << std::endl;

    std::vector<char>::iterator clustersEnd = bclData.cluster(0);
    unsigned clustersLoaded = fastqLoader_.loadClusters(
        bclData.getClusterCount(), flowcell.getReadMetadataList(), clustersEnd);
    ISAAC_ASSERT_MSG(clustersToLoad == clustersLoaded, "Loaded mismatching number of clusters: " << clustersLoaded << " expected: " << tileMetadata);

    bclData.pf().clear();
    for(unsigned cluster = 0; clustersLoaded > cluster; ++cluster)
    {
        bclData.pf().push_back(true);
    }


    ISAAC_THREAD_CERR << "Loading Fastq data done. Loaded " << clustersLoaded << " clusters for " << tileMetadata << std::endl;
}


template class FastqSeedSource<oligo::ShortKmerType>;
template class FastqSeedSource<oligo::KmerType>;
template class FastqSeedSource<oligo::LongKmerType>;

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
