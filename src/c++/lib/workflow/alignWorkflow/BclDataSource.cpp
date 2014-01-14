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
BclSeedSource<KmerT>::BclSeedSource(
    const bool ignoreMissingBcls,
    const unsigned inputLoadersMax,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
    const flowcell::Layout &bclFlowcellLayout,
    common::ThreadVector &threads) :
        ignoreMissingBcls_(ignoreMissingBcls),
        inputLoadersMax_(inputLoadersMax),
        barcodeMetadataList_(barcodeMetadataList),
        bclFlowcellLayout_(bclFlowcellLayout),
        sortedReferenceMetadataList_(sortedReferenceMetadataList),
        flowcellTiles_(getTiles(bclFlowcellLayout)),
        maxTileClusterCount_(std::max_element(flowcellTiles_.begin(), flowcellTiles_.end(),
                                              boost::bind(&flowcell::TileMetadata::getClusterCount, _1)<
                                              boost::bind(&flowcell::TileMetadata::getClusterCount, _2))->getClusterCount()),
        undiscoveredTiles_(flowcellTiles_.begin()),
        threadBclReaders_(inputLoadersMax_,
            rta::BclReader(
                ignoreMissingBcls_,
                maxTileClusterCount_)),
        threads_(threads),
        longestBclPathLength_(
            bclFlowcellLayout_.getLongestAttribute<flowcell::Layout::Bcl, flowcell::BclFilePathAttributeTag>().string().size()),
        barcodeLoader_(threads_, inputLoadersMax_, flowcellTiles_, bclFlowcellLayout_, longestBclPathLength_, threadBclReaders_)
{
    while(threadBclMappers_.size() < inputLoadersMax_)
    {
        threadBclMappers_.push_back(
            new rta::SingleCycleBclMapper<rta::BclReader>(
                maxTileClusterCount_, longestBclPathLength_,
                flowcell::Layout::Bcl == bclFlowcellLayout_.getFormat(),
                threadBclReaders_.at(threadBclMappers_.size())));
    }
}

// TileSource implementation
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

// BarcodeSource implementation
template <typename KmerT>
void BclSeedSource<KmerT>::loadBarcodes(
    const unsigned unknownBarcodeIndex,
    const flowcell::TileMetadataList &tiles,
    std::vector<demultiplexing::Barcode> &barcodes)
{
    barcodeLoader_.loadBarcodes(unknownBarcodeIndex, tiles, barcodes);
}

// SeedSource implementation
template <typename KmerT>
void BclSeedSource<KmerT>::initBuffers(
    flowcell::TileMetadataList &unprocessedTiles,
    const alignment::SeedMetadataList &seedMetadataList)
{
    seedLoader_.reset(new alignment::ParallelSeedLoader<rta::BclReader, KmerT>(
        ignoreMissingBcls_, threads_, threadBclMappers_,
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
            flowcellLayout.getLaneTileCycleAttribute<flowcell::Layout::Bcl, flowcell::BclFilePathAttributeTag>(
                lane, tile, flowcellLayout.getDataCycles().at(0), bclFilePath);
            if (boost::filesystem::exists(bclFilePath))
            {
                const unsigned int clusterCount = rta::BclMapper::getClusterCount(bclFilePath);
                const flowcell::TileMetadata tileMetadata(
                    flowcellId, flowcellLayout.getIndex(),
                    tile, lane,
                    clusterCount,
                    tileMetadataList.size());
                tileMetadataList.push_back(tileMetadata);
            }
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
    const unsigned inputLoadersMax,
    const bool extractClusterXy):
    flowcellLayoutList_(flowcellLayoutList),
    bclLoadThreads_(bclLoadThreads),
    filterFilePath_(flowcell::getLongestAttribute<flowcell::Layout::Bcl, flowcell::FiltersFilePathAttributeTag>(flowcellLayoutList_)),
    positionsFilePath_(flowcell::getLongestAttribute<flowcell::Layout::Bcl, flowcell::PositionsFilePathAttributeTag>(flowcellLayoutList_)),
    threadReaders_(bclLoadThreads_.size(), rta::BclReader(ignoreMissingBcls, flowcell::getMaxTileClusters(tileMetadataList))),
    bclMapper_(ignoreMissingBcls,
               flowcell::getMaxTotalReadLength(flowcellLayoutList_) + flowcell::getMaxBarcodeLength(flowcellLayoutList_),
               bclLoadThreads_, threadReaders_,
               inputLoadersMax, flowcell::getMaxTileClusters(tileMetadataList),
               flowcell::getLongestAttribute<flowcell::Layout::Bcl, flowcell::BclFilePathAttributeTag>(flowcellLayoutList_).string().size()),
    filtersMapper_(ignoreMissingFilters),
    clocsMapper_(),
    locsMapper_()
{
    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions before filtersMapper_.reserveBuffer ")
    filtersMapper_.reserveBuffers(filterFilePath_.string().size(), flowcell::getMaxTileClusters(tileMetadataList));
    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions after filtersMapper_.reserveBuffer ")

    if (extractClusterXy)
    {
        clocsMapper_.reserveBuffers(positionsFilePath_.string().size(), flowcell::getMaxTileClusters(tileMetadataList));
        locsMapper_.reserveBuffers(positionsFilePath_.string().size(), flowcell::getMaxTileClusters(tileMetadataList));
        ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions after filtersMapper_.reserveBuffer ")
    }
}

void BclBaseCallsSource::loadClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData)
{
    ISAAC_THREAD_CERR << "Loading Bcl data for " << tileMetadata << std::endl;

    const flowcell::Layout &flowcell = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex());
    bclMapper_.mapTile(flowcell, tileMetadata);
    ISAAC_THREAD_CERR << "Loading Bcl data done for " << tileMetadata << std::endl;

    ISAAC_THREAD_CERR << "Loading Filter data for " << tileMetadata << std::endl;
    flowcell.getLaneTileAttribute<flowcell::Layout::Bcl, flowcell::FiltersFilePathAttributeTag>(
        tileMetadata.getLane(), tileMetadata.getTile(), filterFilePath_);
    filtersMapper_.mapTile(filterFilePath_, tileMetadata.getClusterCount());
    ISAAC_THREAD_CERR << "Loading Filter data done for " << tileMetadata << std::endl;

    bool boolUseLocsPositions = false;
    if (bclData.storeXy())
    {
        ISAAC_THREAD_CERR << "Loading Positions data for " << tileMetadata << std::endl;
        flowcell.getLaneTileAttribute<flowcell::Layout::Bcl, flowcell::PositionsFilePathAttributeTag>(
            tileMetadata.getLane(), tileMetadata.getTile(), positionsFilePath_);
        if (flowcell::isClocsPath(positionsFilePath_))
        {
            clocsMapper_.mapTile(positionsFilePath_, tileMetadata.getClusterCount());
        }
        else
        {
            locsMapper_.mapTile(positionsFilePath_, tileMetadata.getClusterCount());
            boolUseLocsPositions = true;
        }
        ISAAC_THREAD_CERR << "Loading Positions data done for " << tileMetadata << std::endl;
    }

    // bclToClusters mainly does transposition of bcl cycles to clusters which is a non-io operation.
    // However, the amount of CPU required is relatively low, and occurs on a single thread.
    // Avoid locking all the cores for the duration of this...
    // Also, bclMapper_ and filtersMapper_ are shared between the threads at the moment.
    bclToClusters(tileMetadata, bclData, boolUseLocsPositions);
}


void BclBaseCallsSource::bclToClusters(
    const flowcell::TileMetadata &tileMetadata,
    alignment::BclClusters &bclData,
    const bool useLocsPositions) const
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
        if (!useLocsPositions)
        {
            clocsMapper_.getPositions(std::back_inserter(bclData.xy()));
        }
        else
        {
            locsMapper_.getPositions(std::back_inserter(bclData.xy()));
        }
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
