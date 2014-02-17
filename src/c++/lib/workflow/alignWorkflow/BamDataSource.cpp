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
 ** \file BamDataSource.cpp
 **
 ** \brief see BamDataSource.hh
 **
 ** \author Roman Petrovski
 **/

#include "oligo/Nucleotides.hh"
#include "workflow/alignWorkflow/BamDataSource.hh"
#include "workflow/alignWorkflow/bamDataSource/SingleEndClusterExtractor.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{


/**
 * \param clusterCount  Maximum number of clusters to load
 * \param clusterIt     Insert iterator for the buffer that is sufficient to load the clusterCount
 *                      clusters of clusterLength
 * \param pfIt          Insert iterator for the buffer that is sufficient to load the clusterCount
 *                      clusters of pf flags
 *
 * \return Actual number of loaded clusters
 */
template <typename ClusterInsertIt, typename PfInserIt>
unsigned BamClusterLoader::loadPairedReads(
    unsigned clusterCount,
    const flowcell::ReadMetadataList &readMetadataList,
    ClusterInsertIt &clusterIt, PfInserIt &pfIt)
{
    const unsigned requestedClusterCount = clusterCount;

    while (clusterCount)
    {
        if(clusterExtractor_.extractingUnpaired())
        {
            ISAAC_THREAD_CERR << "extracting unpaired " << std::endl;
            clusterCount = clusterExtractor_.extractUnpaired(
                readMetadataList.at(0).getLength(), readMetadataList.at(1).getLength(), clusterCount,
                clusterIt, pfIt);
            // Either no room in result buffer or no more data available.
            break;
        }
        else
        {
            if (!clusterExtractor_.isEmpty())
            {
                ISAAC_THREAD_CERR << "resuming from " << clusterExtractor_.size() << " pending elements" << std::endl;
                clusterCount = clusterExtractor_.extractPairedReads(clusterCount, clusterIt, pfIt, readMetadataList);

                if (!clusterCount)
                {
                    // there is no room to extract any more pending items. Resume on next call.
                    return requestedClusterCount;
                }
            }

            bamLoader_.load
            (
                boost::make_tuple(
                    boost::bind(&bamDataSource::PairedEndClusterExtractor::append<ClusterInsertIt, PfInserIt>,
                                &clusterExtractor_, _1, _2, boost::ref(clusterCount),
                                boost::ref(readMetadataList), boost::ref(clusterIt), boost::ref(pfIt)),
                    boost::bind(&bamDataSource::PairedEndClusterExtractor::removeOld,
                                &clusterExtractor_, _1, _2, boost::ref(readMetadataList)))
            );


            if (clusterCount)
            {
                // we've ran out of data in the bam file. See if unpaired items can be paired
                clusterExtractor_.startExtractingUnpaired();
            }
        }
    }

    ISAAC_THREAD_CERR << "loadPairedReads done clusterCount:" << clusterCount << std::endl;
    return requestedClusterCount - clusterCount;
}

void BamClusterLoader::open(
    const std::string &flowcellId,
    const boost::filesystem::path &bamPath)
{
    if (flowcellId_ != flowcellId)
    {
        clusterExtractor_.open(flowcellId);
        bamLoader_.open(bamPath);
        flowcellId_ = flowcellId;
    }
    else
    {
        ISAAC_THREAD_CERR << "Keeping bam stream open for flowcellId " << flowcellId_ << " " << bamPath << std::endl;
    }
}

template <typename InsertIt, typename PfInserIt>
unsigned BamClusterLoader::loadSingleReads(
    unsigned clusterCount, const flowcell::ReadMetadataList &readMetadataList,
    InsertIt &clusterIt, PfInserIt &pfIt)
{
    const unsigned requestedClusterCount = clusterCount;
    bamDataSource::SingleEndClusterExtractor extractor;
    bamLoader_.load
    (
        boost::make_tuple(
            boost::bind(&bamDataSource::SingleEndClusterExtractor::extractSingleRead<InsertIt, PfInserIt>,
                &extractor, _1, boost::ref(clusterCount),
                boost::ref(readMetadataList), boost::ref(clusterIt), boost::ref(pfIt)),
            boost::bind(&bamDataSource::SingleEndClusterExtractor::nothing))
    );

    ISAAC_THREAD_CERR << "loadSingleReads done clusterCount:" << clusterCount << " bgzfReader_.isEof():" << std::endl;

    return requestedClusterCount - clusterCount;}

template <typename ClusterInsertIt, typename PfInserIt>
unsigned BamClusterLoader::loadClusters(
    unsigned clusterCount, const flowcell::ReadMetadataList &readMetadataList,
    ClusterInsertIt &clusterIt, PfInserIt &pfIt)
{
    return 2 == readMetadataList.size() ?
        loadPairedReads(clusterCount, readMetadataList, clusterIt, pfIt) :
        loadSingleReads(clusterCount, readMetadataList, clusterIt, pfIt);
}

template <typename KmerT>
unsigned BamSeedSource<KmerT>::determineMemoryCapacity(
    const unsigned long availableMemory,
    const unsigned tileClustersMax,
    const unsigned clusterLength)
{
    ISAAC_THREAD_CERR << "Determining memory capacity for Bam data" << std::endl;
    // Allocate as much memory as possible in tileClustersMax_ decrements
    alignment::BclClusters testClusters(clusterLength);
    // Assume seeds will equate at most to the cluster length at a time. This
    // means the RAM available for clusters will be third (seeds are made of kmer and metadata)
    // the ram successfully allocated
    size_t testCount = availableMemory / clusterLength;
    while (testCount)
    {
        try
        {
            ISAAC_THREAD_CERR << "Determining memory capacity trying. " << testCount << std::endl;
            testClusters.reserveClusters(testCount, false);
            break;
        }
        catch (std::bad_alloc &e)
        {
            if (tileClustersMax >= testCount)
            {
                BOOST_THROW_EXCEPTION(common::MemoryException((boost::format(
                    "Insufficient memory to load seeds even for %u bam clusters") % testCount).str()));
            }
            // reset errno, to prevent misleading error messages when failing code does not set errno
            errno = 0;
            testCount -= tileClustersMax;
        }
    }

    ISAAC_THREAD_CERR << "Determining memory capacity for Bam data done. " <<
        testCount << " clusters of length " << clusterLength <<
            " will fit" << std::endl;
    return testCount;
}


inline std::size_t getBamFileSize(const flowcell::Layout &flowcell)
{
    return common::getFileSize(flowcell.getAttribute<flowcell::Layout::Bam, flowcell::BamFilePathAttributeTag>().c_str());
}

template <typename KmerT>
BamSeedSource<KmerT>::BamSeedSource(
    const boost::filesystem::path &tempDirectoryPath,
    const unsigned long availableMemory,
    const bool cleanupIntermediary,
    const unsigned coresMax,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
    const flowcell::Layout &bamFlowcellLayout,
    common::ThreadVector &threads) :
        bamFlowcellLayout_(bamFlowcellLayout),
        tileClustersMax_(40000000 / bamFlowcellLayout_.getSeedMetadataList().size()),
        coresMax_(coresMax),
        barcodeMetadataList_(barcodeMetadataList),
        sortedReferenceMetadataList_(sortedReferenceMetadataList),
        clusterLength_(flowcell::getTotalReadLength(bamFlowcellLayout_.getReadMetadataList())),
        clustersAtATimeMax_(determineMemoryCapacity(availableMemory, tileClustersMax_, clusterLength_)),
        clusters_(clusterLength_),
        currentTile_(1),
        threads_(threads),
        bamClusterLoader_(
            cleanupIntermediary, 0, threads, coresMax, tempDirectoryPath,
            getBamFileSize(bamFlowcellLayout_), bamFlowcellLayout.getFlowcellId().length(),
            flowcell::getTotalReadLength(bamFlowcellLayout.getReadMetadataList()))
{
}

template <typename T>
struct VoidInsertIterator
{
    VoidInsertIterator operator ++() const
    {
        return VoidInsertIterator();
    }
    VoidInsertIterator operator ++(int) const
    {
        return VoidInsertIterator();
    }

    VoidInsertIterator operator =(const T &that) const
    {
        return VoidInsertIterator();
    }

    VoidInsertIterator operator *() const
    {
        return VoidInsertIterator();
    }
};

template <typename KmerT>
flowcell::TileMetadataList BamSeedSource<KmerT>::discoverTiles()
{
    loadedTiles_.clear();

    // Allocate one third the available RAM for clusters as for each cluster kmer we need a kmer+metadata for the seed
    const unsigned clustersToLoad = clustersAtATimeMax_ / 3;
    clusters_.reset(clusterLength_, clustersToLoad);
    // load clusters, return tile breakdown based on tileClustersMax_

    const boost::filesystem::path bamPath = bamFlowcellLayout_.getAttribute<flowcell::Layout::Bam, flowcell::BamFilePathAttributeTag>();
    // this will keep the current files open if the paths don't change
    bamClusterLoader_.open(bamFlowcellLayout_.getFlowcellId(), bamPath);

    std::vector<char>::iterator clustersEnd = clusters_.cluster(0);
    VoidInsertIterator<bool> dummy;
    unsigned clustersLoaded = bamClusterLoader_.loadClusters(
        clustersToLoad, bamFlowcellLayout_.getReadMetadataList(), clustersEnd, dummy);
    ISAAC_THREAD_CERR << "Loaded  " << clustersLoaded << " clusters of length " << clusterLength_ << std::endl;
    clusters_.reset(clusterLength_, clustersLoaded);
    if (clustersLoaded)
    {
        if (clustersLoaded < clustersToLoad / 2)
        {
            //allocated too much memory for bcl data. Might pay off to give some away especially if seeds are placed very densely for high sensitivity
            clusters_.reduceWastedMemory();
        }
        const std::string &flowcellId = bamFlowcellLayout_.getFlowcellId();
        std::vector<char>::iterator tileFirstCluster = clusters_.cluster(0);
        while (true)
        {
            const unsigned clusterCount = std::min(clustersLoaded, tileClustersMax_);
            const flowcell::TileMetadata tileMetadata(
                flowcellId, bamFlowcellLayout_.getIndex(),
                currentTile_++, 1,
                clusterCount,
                loadedTiles_.size());
            loadedTiles_.push_back(tileMetadata);
            ISAAC_THREAD_CERR << "Generated bam tile: " << oligo::bclToString(reinterpret_cast<const unsigned char *>(&*tileFirstCluster), 100) << "(" <<
                oligo::bclToRString(reinterpret_cast<const unsigned char *>(&*tileFirstCluster), 100) << ")" << tileMetadata << std::endl;
            tileFirstCluster += clusterCount * flowcell::getTotalReadLength(bamFlowcellLayout_.getReadMetadataList());

            if (clustersLoaded < tileClustersMax_)
            {
                break;
            }
            clustersLoaded -= tileClustersMax_;
        }
    }

    return loadedTiles_;
}

template <typename KmerT>
void BamSeedSource<KmerT>::initBuffers(
    flowcell::TileMetadataList &unprocessedTiles,
    const alignment::SeedMetadataList &seedMetadataList)
{
    seedGenerator_.reset(new alignment::ClusterSeedGenerator<KmerT>(
        threads_,
        coresMax_, barcodeMetadataList_,
        bamFlowcellLayout_,
        seedMetadataList,
        sortedReferenceMetadataList_, unprocessedTiles,
        clusters_,
        loadedTiles_));
}

template <typename KmerT>
void BamSeedSource<KmerT>::generateSeeds(
    const flowcell::TileMetadataList &tiles,
    const alignment::matchFinder::TileClusterInfo &tileClusterBarcode,
    std::vector<SeedT> &seeds,
    common::ScoopedMallocBlock  &mallocBlock)
{
    seedGenerator_->generateSeeds(tiles, tileClusterBarcode, seeds, mallocBlock);
}

template <typename KmerT>
const std::vector<typename BamSeedSource<KmerT>::SeedIterator> &BamSeedSource<KmerT>::getReferenceSeedBounds() const
{
    return seedGenerator_->getReferenceSeedBounds();
}

/////////////// BamBaseCallsSource implementation
inline const boost::filesystem::path getLongestBamFilePath(const flowcell::FlowcellLayoutList &flowcellLayoutList)
{
    boost::filesystem::path ret;
    BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
    {
        if (flowcell::Layout::Bam == flowcell.getFormat())
        {
            const boost::filesystem::path bamPath = flowcell.getAttribute<flowcell::Layout::Bam, flowcell::BamFilePathAttributeTag>();
            if (bamPath.string().length() > ret.string().length())
            {
                ret = bamPath;
            }
        }
    }
    return ret;
}

inline std::size_t getBiggestBamFileSize(const flowcell::FlowcellLayoutList &flowcellLayoutList)
{
    std::size_t ret = 0;
    BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
    {
        if (flowcell::Layout::Bam == flowcell.getFormat())
        {
            const std::size_t bamSize = getBamFileSize(flowcell);
            ret = std::max(ret, bamSize);
        }
    }
    return ret;
}

inline std::size_t getLongestFlowcellId(const flowcell::FlowcellLayoutList &flowcellLayoutList)
{
    return std::max_element(flowcellLayoutList.begin(), flowcellLayoutList.end(),
        boost::bind(&std::string::size, boost::bind(&flowcell::Layout::getFlowcellId, _1)) <
        boost::bind(&std::string::size, boost::bind(&flowcell::Layout::getFlowcellId, _2))
    )->getFlowcellId().length();
}


BamBaseCallsSource::BamBaseCallsSource(
    const boost::filesystem::path &tempDirectoryPath,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const flowcell::TileMetadataList &tileMetadataList,
    const bool cleanupIntermediary,
    common::ThreadVector &threads,
    const unsigned inputLoadersMax):
    flowcellLayoutList_(flowcellLayoutList),
    bamClusterLoader_(
        cleanupIntermediary,
        getLongestBamFilePath(flowcellLayoutList_).string().length(),
        threads, inputLoadersMax, tempDirectoryPath,
        getBiggestBamFileSize(flowcellLayoutList),
        getLongestFlowcellId(flowcellLayoutList_),
        flowcell::getMinTotalReadLength(flowcellLayoutList_)),
    // reserve memory for auxiliary structures needed for processing
    bamFilePath_(getLongestBamFilePath(flowcellLayoutList_).c_str())
{
}

void BamBaseCallsSource::loadClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData)
{
    ISAAC_THREAD_CERR << "Loading Bam data for " << tileMetadata << std::endl;
    const flowcell::Layout &flowcell = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex());
    ISAAC_ASSERT_MSG(1 == tileMetadata.getLane(), "Bam files are expected to have only one 'lane'");

    flowcell.getAttribute<flowcell::Layout::Bam, flowcell::BamFilePathAttributeTag>(bamFilePath_);

    bamClusterLoader_.open(flowcell.getFlowcellId(), bamFilePath_);

    const unsigned clustersToLoad = tileMetadata.getClusterCount();
    ISAAC_THREAD_CERR << "Resetting Bam data for " << clustersToLoad << " clusters" << std::endl;
    bclData.reset(flowcell::getTotalReadLength(flowcell.getReadMetadataList()), clustersToLoad);
    ISAAC_THREAD_CERR << "Resetting Bam data done for " << bclData.getClusterCount() << " clusters" << std::endl;

    std::vector<char>::iterator clusterIt = bclData.cluster(0);
    bclData.pf().clear();
    std::back_insert_iterator<std::vector<bool> > pfIt(bclData.pf());
    unsigned clustersLoaded = bamClusterLoader_.loadClusters(
        bclData.getClusterCount(), flowcell.getReadMetadataList(), clusterIt, pfIt);
    ISAAC_ASSERT_MSG(clustersToLoad == clustersLoaded, "Loaded mismatching number of clusters: " << clustersLoaded << " expected: " << tileMetadata);

//    ISAAC_THREAD_CERR << "Loaded bam tile: " << oligo::bclToString(reinterpret_cast<const unsigned char *>(&*bclData.cluster(0)), 100) << "(" <<
//        oligo::bclToRString(reinterpret_cast<const unsigned char *>(&*bclData.cluster(0)), 100) << ")" << tileMetadata << std::endl;

    ISAAC_THREAD_CERR << "Loading Bam data done. Loaded " << clustersLoaded << " clusters for " << tileMetadata << std::endl;
}


template class BamSeedSource<oligo::ShortKmerType>;
template class BamSeedSource<oligo::KmerType>;
template class BamSeedSource<oligo::LongKmerType>;

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
