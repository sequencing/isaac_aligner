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
 ** \file BclDataSource.cpp
 **
 ** \brief see BclDataSource.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>

#include "workflow/alignWorkflow/BclDataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{


template <typename KmerT>
flowcell::TileMetadataList BclSeedSource<KmerT>::discoverTiles()
{
    flowcell::TileMetadataList ret;
    while (flowcellTiles_.end() != undiscoveredTiles_)
    {
        ret.push_back(*undiscoveredTiles_);
        if (ret.front().getLane() != ret.back().getLane())
        {
            ret.pop_back();
            break;
        }
        ++undiscoveredTiles_;
    }
    return ret;
}

template <typename KmerT>
void BclSeedSource<KmerT>::initBuffers(
    flowcell::TileMetadataList &unprocessedTiles,
    const alignment::SeedMetadataList &seedMetadataList,
    common::ThreadVector &threads)
{
    seedLoader_.reset(new alignment::ParallelSeedLoader<KmerT>(
        ignoreMissingBcls_, threads,
        inputLoadersMax_, barcodeMetadataList_,
        bclFlowcellLayout_,
        seedMetadataList,
        sortedReferenceMetadataList_, unprocessedTiles));
}

template <typename KmerT>
void BclSeedSource<KmerT>::generateSeeds(
    const flowcell::TileMetadataList &tiles,
    const alignment::matchFinder::TileClusterInfo &tileClusterBarcode,
    std::vector<SeedT> &seeds,
    common::ScoopedMallocBlock  &mallocBlock)
{
    seedLoader_->loadSeeds(tiles, tileClusterBarcode, seeds, mallocBlock);
}

template <typename KmerT>
const std::vector<typename BclSeedSource<KmerT>::SeedIterator> &BclSeedSource<KmerT>::getReferenceSeedBounds() const
{
    return seedLoader_->getReferenceSeedBounds();
}

template <typename KmerT>
flowcell::TileMetadataList BclSeedSource<KmerT>::getTiles(const flowcell::Layout &flowcellLayout) const
{
    flowcell::TileMetadataList tileMetadataList;

    const std::string &flowcellId = flowcellLayout.getFlowcellId();
    BOOST_FOREACH(const unsigned int lane, flowcellLayout.getLaneIds())
    {
        const std::vector<unsigned int> tileList = flowcellLayout.getTileIds(lane);
        BOOST_FOREACH(const unsigned int tile, tileList)
        {
            boost::filesystem::path bclFilePath;
            flowcellLayout.getBclFilePath(tile, lane, flowcellLayout.getAllCycleNumbers().at(0), bclFilePath);
            const unsigned int clusterCount = io::BclMapper::getClusterCount(bclFilePath);
            const flowcell::TileMetadata tileMetadata(
                flowcellId, flowcellLayout.getIndex(),
                tile, lane,
                clusterCount,
                tileMetadataList.size());
            tileMetadataList.push_back(tileMetadata);
        }
    }

    return tileMetadataList;
}

/////////////// BclBaseCallsSource Implementation

BclBaseCallsSource::BclBaseCallsSource(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const flowcell::TileMetadataList &tileMetadataList,
    const bool ignoreMissingBcls,
    const bool ignoreMissingFilters,
    common::ThreadVector &bclLoadThreads,
    const unsigned inputLoadersMax):
    flowcellLayoutList_(flowcellLayoutList),
    bclLoadThreads_(bclLoadThreads),
    bclFilePaths_(flowcell::getMaxTotalReadLength(flowcellLayoutList_) + flowcell::getMaxBarcodeLength(flowcellLayoutList_)),
    filterFilePath_(flowcell::Layout::getLongestFilterFilePath(flowcellLayoutList_)),
    positionsFilePath_(flowcell::Layout::getLongestPositionsFilePath(flowcellLayoutList_)),
    bclMapper_(ignoreMissingBcls, flowcell::getMaxTotalReadLength(flowcellLayoutList_) +
           flowcell::getMaxBarcodeLength(flowcellLayoutList_),
           bclLoadThreads_, inputLoadersMax, flowcell::getMaxTileClusters(tileMetadataList)),
    filtersMapper_(ignoreMissingFilters),
    positionsMapper_()
{
    // reserve memory for auxiliary structures needed for bcl processing
    const boost::filesystem::path longestBaseCallsPath =  flowcell::getLongestBaseCallsPath(flowcellLayoutList_);
    const bool compressedFound =
        flowcellLayoutList_.end() != std::find_if(flowcellLayoutList_.begin(), flowcellLayoutList_.end(),
                                                boost::bind(&flowcell::Layout::getFormat, _1) ==
                                                    flowcell::Layout::BclGz);

    const unsigned highestTileNumber = std::max_element(tileMetadataList.begin(), tileMetadataList.end(),
                                                        boost::bind(&flowcell::TileMetadata::getTile, _1)<
                                                            boost::bind(&flowcell::TileMetadata::getTile, _2))->getTile();
    const unsigned insanelyHighCycleNumber = 9999;
    boost::filesystem::path longestBclFilePath;
    flowcell::Layout::getBclFilePath(highestTileNumber, 1, longestBaseCallsPath, insanelyHighCycleNumber,
                                     true, longestBclFilePath);

    bclMapper_.reserveBuffers(longestBclFilePath.string().size(), compressedFound);

    BOOST_FOREACH(boost::filesystem::path &p, bclFilePaths_)
    {
        // this has to be done separately for each path or else they all share one buffer
        p = longestBclFilePath.c_str();
    }

    filtersMapper_.reservePathBuffers(filterFilePath_.string().size());
    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions before filtersMapper_.reserveBuffer ")
    filtersMapper_.reserveBuffer(flowcell::getMaxTileClusters(tileMetadataList));
    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions after filtersMapper_.reserveBuffer ")

    positionsMapper_.reservePathBuffers(positionsFilePath_.string().size());
    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions before filtersMapper_.reserveBuffer ")
    positionsMapper_.reserveBuffer(flowcell::getMaxTileClusters(tileMetadataList));
    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions after filtersMapper_.reserveBuffer ")
}

void BclBaseCallsSource::loadClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData)
{
    ISAAC_THREAD_CERR << "Loading Bcl data for " << tileMetadata << std::endl;

    const flowcell::Layout &flowcell = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex());
    const std::vector<unsigned> &cycleList = flowcell.getAllCycleNumbers();
    ISAAC_ASSERT_MSG(bclFilePaths_.size() >= cycleList.size() + flowcell.getBarcodeCycles().size(), "tiles expected to be ordered so that number of cycles goes down")
    bclFilePaths_.resize(cycleList.size() + flowcell.getBarcodeCycles().size());
    unsigned pos = 0;

    // Add barcodes first
    BOOST_FOREACH(const unsigned int cycle, flowcell.getBarcodeCycles())
    {
        flowcellLayoutList_.at(tileMetadata.getFlowcellIndex()).getBclFilePath(
            tileMetadata.getTile(), tileMetadata.getLane(),
            cycle, bclFilePaths_[pos++]);
    }

    // Add reads (non-barcodes) second
    BOOST_FOREACH(const unsigned int cycle, cycleList)
    {
        flowcellLayoutList_.at(tileMetadata.getFlowcellIndex()).getBclFilePath(
            tileMetadata.getTile(), tileMetadata.getLane(),
            cycle, bclFilePaths_[pos++]);
    }

    bclMapper_.mapTile(bclFilePaths_, tileMetadata.getClusterCount());
    ISAAC_THREAD_CERR << "Loading Bcl data done for " << tileMetadata << std::endl;

    ISAAC_THREAD_CERR << "Loading Filter data for " << tileMetadata << std::endl;
    flowcell.getFiltersFilePath(tileMetadata.getTile(), tileMetadata.getLane(), filterFilePath_);
    filtersMapper_.mapTile(filterFilePath_, tileMetadata.getClusterCount());
    ISAAC_THREAD_CERR << "Loading Filter data done for " << tileMetadata << std::endl;

    if (bclData.storeXy())
    {
        ISAAC_THREAD_CERR << "Loading Positions data for " << tileMetadata << std::endl;
        flowcell.getPositionsFilePath(tileMetadata.getTile(), tileMetadata.getLane(), positionsFilePath_);
        positionsMapper_.mapTile(positionsFilePath_, tileMetadata.getClusterCount());
        ISAAC_THREAD_CERR << "Loading Positions data done for " << tileMetadata << std::endl;
    }

    // bclToClusters mainly does transposition of bcl cycles to clusters which is a non-io operation.
    // However, the amount of CPU required is relatively low, and occurs on a single thread.
    // Avoid locking all the cores for the duration of this...
    // Also, bclMapper_ and filtersMapper_ are shared between the threads at the moment.
    bclToClusters(tileMetadata, bclData);
}


void BclBaseCallsSource::bclToClusters(
    const flowcell::TileMetadata &tileMetadata,
    alignment::BclClusters &bclData) const
{

    ISAAC_THREAD_CERR << "Resetting Bcl data for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
    bclData.reset(bclMapper_.getCyclesCount(), tileMetadata.getClusterCount());
    ISAAC_THREAD_CERR << "Resetting Bcl data done for " << bclData.getClusterCount() << " bcl clusters" << std::endl;


    ISAAC_THREAD_CERR << "Transposing Bcl data for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
    const clock_t startTranspose = clock();
    bclMapper_.transpose(bclData.cluster(0));
    ISAAC_THREAD_CERR << "Transposing Bcl data done for " << bclData.getClusterCount() << " bcl clusters in " << (clock() - startTranspose) / 1000 << "ms" << std::endl;

    ISAAC_THREAD_CERR << "Extracting Pf values for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
    bclData.pf().clear();
    filtersMapper_.getPf(std::back_inserter(bclData.pf()));
    assert(bclData.pf().size() == tileMetadata.getClusterCount());
    ISAAC_THREAD_CERR << "Extracting Pf values done for " << bclData.getClusterCount() << " bcl clusters" << std::endl;

    if (bclData.storeXy())
    {
        ISAAC_THREAD_CERR << "Extracting Positions values for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
        bclData.xy().clear();
        positionsMapper_.getPositions(std::back_inserter(bclData.xy()));
        assert(bclData.xy().size() == tileMetadata.getClusterCount());
        ISAAC_THREAD_CERR << "Extracting Positions values done for " << bclData.getClusterCount() << " bcl clusters" << std::endl;
    }

}

template class BclSeedSource<oligo::ShortKmerType>;
template class BclSeedSource<oligo::KmerType>;
template class BclSeedSource<oligo::LongKmerType>;

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
