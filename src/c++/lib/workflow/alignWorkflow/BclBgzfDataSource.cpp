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
 ** \file BclDataSource.cpp
 **
 ** \brief see BclDataSource.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>

#include "workflow/alignWorkflow/BclBgzfDataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

template <typename KmerT>
BclBgzfSeedSource<KmerT>::BclBgzfSeedSource(
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
        threads_(threads),
        sortedReferenceMetadataList_(sortedReferenceMetadataList),
        tileBciIndexMap_(),
        flowcellTiles_(getTiles(bclFlowcellLayout_, tileBciIndexMap_)),
        maxTileClusterCount_(std::max_element(flowcellTiles_.begin(), flowcellTiles_.end(),
                                              boost::bind(&flowcell::TileMetadata::getClusterCount, _1)<
                                              boost::bind(&flowcell::TileMetadata::getClusterCount, _2))->getClusterCount()),
        undiscoveredTiles_(flowcellTiles_.begin()),
        threadBclReaders_(inputLoadersMax_,
            rta::BclBgzfTileReader(ignoreMissingBcls_, maxTileClusterCount_, tileBciIndexMap_, cycleBciMappers_)),
        barcodeLoader_(
            threads_, inputLoadersMax_, flowcellTiles_,
            bclFlowcellLayout_,
            bclFlowcellLayout_.getLongestAttribute<flowcell::Layout::BclBgzf, flowcell::BclFilePathAttributeTag>().string().size(),
            threadBclReaders_)
{
    const unsigned longestBclPath =
        bclFlowcellLayout_.getLongestAttribute<flowcell::Layout::BclBgzf, flowcell::BclFilePathAttributeTag>().string().size();

    while(threadBclMappers_.size() < inputLoadersMax_)
    {
        threadBclMappers_.push_back(
            new rta::SingleCycleBclMapper<rta::BclBgzfTileReader>(
                maxTileClusterCount_, longestBclPath,
                true,
                threadBclReaders_.at(threadBclMappers_.size())));
    }
}

// TileSource implementation
template <typename KmerT>
flowcell::TileMetadataList BclBgzfSeedSource<KmerT>::discoverTiles()
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

inline unsigned getLaneNumber(const flowcell::TileMetadataList& tiles)
{
    ISAAC_ASSERT_MSG(
        tiles.end()
            == std::find_if(
                tiles.begin(),
                tiles.end(),
                boost::bind(&flowcell::TileMetadata::getLane, _1)
                    != tiles.front().getLane()),
        "Expected all tiles to belong to the same lane");
    return tiles.front().getLane();
}

// BarcodeSource implementation
template <typename KmerT>
void BclBgzfSeedSource<KmerT>::loadBarcodes(
    const unsigned unknownBarcodeIndex,
    const flowcell::TileMetadataList &tiles,
    std::vector<demultiplexing::Barcode> &barcodes)
{
    const std::vector<unsigned> cycles = bclFlowcellLayout_.getBarcodeCycles();
    initCycleBciMappers(cycles, getLaneNumber(tiles), cycleBciMappers_);

    barcodeLoader_.loadBarcodes(unknownBarcodeIndex, tiles, barcodes);
}

template<typename KmerT>
void BclBgzfSeedSource<KmerT>::initCycleBciMappers(
    const std::vector<unsigned>& cycles,
    const unsigned laneNumber,
    std::vector<rta::CycleBciMapper> &cycleBciMappers) const
{
    cycleBciMappers.resize(cycles.back() + 1, rta::CycleBciMapper(0));
    BOOST_FOREACH(const unsigned cycle, cycles)
    {
        boost::filesystem::path cycleBciFilePath;
        bclFlowcellLayout_.getLaneCycleAttribute<
            flowcell::Layout::BclBgzf,flowcell::BciFilePathAttributeTag>(laneNumber, cycle, cycleBciFilePath);
        cycleBciMappers.at(cycle).mapFile(cycleBciFilePath);
    }
}

// SeedSource implementation
template <typename KmerT>
void BclBgzfSeedSource<KmerT>::initBuffers(
    flowcell::TileMetadataList &unprocessedTiles,
    const alignment::SeedMetadataList &seedMetadataList)
{
    const std::vector<unsigned> cycles = alignment::getAllSeedCycles(
        bclFlowcellLayout_.getReadMetadataList(), seedMetadataList);
    initCycleBciMappers(cycles, getLaneNumber(unprocessedTiles), cycleBciMappers_);

    seedLoader_.reset(new alignment::ParallelSeedLoader<rta::BclBgzfTileReader, KmerT>(
        ignoreMissingBcls_, threads_, threadBclMappers_,
        inputLoadersMax_, barcodeMetadataList_,
        bclFlowcellLayout_,
        seedMetadataList,
        sortedReferenceMetadataList_, unprocessedTiles));
}

template <typename KmerT>
void BclBgzfSeedSource<KmerT>::generateSeeds(
    const flowcell::TileMetadataList &tiles,
    const alignment::matchFinder::TileClusterInfo &tileClusterBarcode,
    std::vector<SeedT> &seeds,
    common::ScoopedMallocBlock  &mallocBlock)
{
    seedLoader_->loadSeeds(tiles, tileClusterBarcode, seeds, mallocBlock);
}

template <typename KmerT>
const std::vector<typename BclBgzfSeedSource<KmerT>::SeedIterator> &BclBgzfSeedSource<KmerT>::getReferenceSeedBounds() const
{
    return seedLoader_->getReferenceSeedBounds();
}

template <typename KmerT>
flowcell::TileMetadataList BclBgzfSeedSource<KmerT>::getTiles(
    const flowcell::Layout &flowcellLayout,
    std::vector<unsigned> &tileBciIndexMap)
{
    flowcell::TileMetadataList tileMetadataList;

    const std::string &flowcellId = flowcellLayout.getFlowcellId();
    BOOST_FOREACH(const unsigned int lane, flowcellLayout.getLaneIds())
    {
        const std::vector<unsigned int> tileList = flowcellLayout.getTileIds(lane);

        boost::filesystem::path laneBciFilePath;
        flowcellLayout.getLaneAttribute<flowcell::Layout::BclBgzf, flowcell::BciFilePathAttributeTag>(lane, laneBciFilePath);
        rta::LaneBciMapper mapper(
            flowcellLayout.getAttribute<flowcell::Layout::BclBgzf, flowcell::TilesPerLaneMaxAttributeTag>());
        mapper.mapFile(laneBciFilePath);

        BOOST_FOREACH(const unsigned int tileNumber, tileList)
        {
            const rta::LaneBciMapper::TileInfo tileInfo = mapper.getTileInfo(tileNumber);
            if (tileInfo.tileClusters_)
            {
                const flowcell::TileMetadata tileMetadata(
                    flowcellId, flowcellLayout.getIndex(),
                    tileNumber, lane,
                    tileInfo.tileClusters_,
                    tileMetadataList.size());
                tileMetadataList.push_back(tileMetadata);
                tileBciIndexMap.resize(std::max<unsigned>(tileBciIndexMap.size(), tileMetadata.getIndex() + 1));
                tileBciIndexMap[tileMetadata.getIndex()] = tileInfo.tileIndex_;
                ISAAC_THREAD_CERR << tileMetadata << std::endl;
            }
        }
    }

    return tileMetadataList;
}

/////////////// BclBgzfBaseCallsSource Implementation

BclBgzfBaseCallsSource::BclBgzfBaseCallsSource(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const flowcell::TileMetadataList &tileMetadataList,
    const bool ignoreMissingBcls,
    const bool ignoreMissingFilters,
    common::ThreadVector &bclLoadThreads,
    const unsigned inputLoadersMax,
    const bool extractClusterXy):
    flowcellLayoutList_(flowcellLayoutList),
    bclLoadThreads_(bclLoadThreads),
    bciFilePath_(flowcell::getLongestAttribute<flowcell::Layout::BclBgzf, flowcell::BciFilePathAttributeTag>(flowcellLayoutList_)),
    filterFilePath_(flowcell::getLongestAttribute<flowcell::Layout::BclBgzf, flowcell::FiltersFilePathAttributeTag>(flowcellLayoutList_)),
    positionsFilePath_(flowcell::getLongestAttribute<flowcell::Layout::BclBgzf, flowcell::PositionsFilePathAttributeTag>(flowcellLayoutList_)),
    cycles_(flowcell::getMaxTotalReadLength(flowcellLayoutList_) + flowcell::getMaxBarcodeLength(flowcellLayoutList_)),
    cycleBciFileBuf_(std::ios_base::binary|std::ios_base::in),
    tileClusterOffsets_(flowcell::getMaxAttribute<flowcell::Layout::BclBgzf, flowcell::TilesPerLaneMaxAttributeTag>(flowcellLayoutList)),
    cycleBciMappers_(
        flowcell::getMaxCycleNumber(flowcellLayoutList) + 1,
        rta::CycleBciMapper(tileClusterOffsets_.size())),
    tileBciIndexMap_(tileMetadataList.size()),
    laneBciMapper_(tileClusterOffsets_.size()),
    threadReaders_(
        bclLoadThreads_.size(),
        rta::BclBgzfTileReader(ignoreMissingBcls,
                               flowcell::getMaxTileClusters(tileMetadataList),
                               tileBciIndexMap_,
                               cycleBciMappers_)),
    bclMapper_(ignoreMissingBcls, cycles_.size(),
           bclLoadThreads_, threadReaders_,
           inputLoadersMax, flowcell::getMaxTileClusters(tileMetadataList),
           flowcell::getLongestAttribute<flowcell::Layout::BclBgzf, flowcell::BclFilePathAttributeTag>(flowcellLayoutList_).string().size()),
    filtersMapper_(ignoreMissingFilters),
    clocsMapper_(),
    locsMapper_(),
    currentFlowcellIndex_(-1U),
    currentLaneNumber_(-1U)
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

void BclBgzfBaseCallsSource::initBciMappers(
    const flowcell::Layout &flowcell,
    const flowcell::TileMetadataList &allTiles,
    const unsigned laneNumber,
    std::vector<rta::CycleBciMapper>  &cycleBciMappers,
    std::vector<unsigned long> &tileClusterOffsets,
    std::vector<unsigned> &tileBciIndexMap)
{
    flowcell.getLaneAttribute<flowcell::Layout::BclBgzf, flowcell::BciFilePathAttributeTag>(laneNumber, bciFilePath_);
    laneBciMapper_.mapFile(bciFilePath_);
    tileBciIndexMap.resize(allTiles.size());
    BOOST_FOREACH(const flowcell::TileMetadata &tileMetadata, allTiles)
    {
        if (tileMetadata.getFlowcellIndex() == flowcell.getIndex() && laneNumber == tileMetadata.getLane())
        {
            const rta::LaneBciMapper::TileInfo tileInfo = laneBciMapper_.getTileInfo(tileMetadata.getTile());
            tileBciIndexMap.at(tileMetadata.getIndex()) = tileInfo.tileIndex_;
        }
    }

    unsigned prev = 0;
    tileClusterOffsets.resize(laneBciMapper_.getTilesCount());
    for (unsigned bciTileIndex = 0; bciTileIndex < tileClusterOffsets.size(); ++bciTileIndex)
    {
        tileClusterOffsets.at(bciTileIndex) = prev;
        prev += laneBciMapper_.getTileClusterCount(bciTileIndex);
    }

    cycles_.clear();
    cycles_.insert(cycles_.end(), flowcell.getBarcodeCycles().begin(), flowcell.getBarcodeCycles().end());
    cycles_.insert(cycles_.end(), flowcell.getDataCycles().begin(), flowcell.getDataCycles().end());

    ISAAC_ASSERT_MSG(cycleBciMappers.size() > *std::max_element(cycles_.begin(), cycles_.end()),
                     "Insufficient number of cycleBciMappers. Expected: >" <<
                     *std::max_element(cycles_.begin(), cycles_.end()) << " have:" << cycleBciMappers.size());
    BOOST_FOREACH(const unsigned cycle, cycles_)
    {
        flowcell.getLaneCycleAttribute<
            flowcell::Layout::BclBgzf,flowcell::BciFilePathAttributeTag>(laneNumber, cycle, bciFilePath_);

        std::istream is(cycleBciFileBuf_.reopen(bciFilePath_.c_str(), io::FileBufWithReopen::SequentialOnce));
        if (!is)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open file: " + bciFilePath_.string() + strerror(errno)));
        }

        cycleBciMappers.at(cycle).mapStream(is, bciFilePath_);
    }
}

void BclBgzfBaseCallsSource::loadClusters(
    const flowcell::TileMetadataList &allTiles,
    const flowcell::TileMetadata &tileMetadata,
    alignment::BclClusters &bclData)
{
    ISAAC_THREAD_CERR << "Loading Bcl data for " << tileMetadata << std::endl;

    const flowcell::Layout &flowcell = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex());

    if (currentFlowcellIndex_ != flowcell.getIndex() || currentLaneNumber_ != tileMetadata.getLane())
    {
        initBciMappers(flowcell, allTiles, tileMetadata.getLane(), cycleBciMappers_, tileClusterOffsets_, tileBciIndexMap_);
        currentFlowcellIndex_ = flowcell.getIndex();
        currentLaneNumber_ = tileMetadata.getLane();
    }

    bclMapper_.mapTile(flowcell, tileMetadata);
    ISAAC_THREAD_CERR << "Loading Bcl data done for " << tileMetadata << std::endl;

    ISAAC_THREAD_CERR << "Loading Filter data for " << tileMetadata << std::endl;
    flowcell.getLaneAttribute<flowcell::Layout::BclBgzf, flowcell::FiltersFilePathAttributeTag>(tileMetadata.getLane(), filterFilePath_);
    filtersMapper_.mapTile(filterFilePath_, tileMetadata.getClusterCount(),
                           tileClusterOffsets_.at(tileBciIndexMap_.at(tileMetadata.getIndex())));
    ISAAC_THREAD_CERR << "Loading Filter data done for " << tileMetadata << std::endl;

    bool boolUseLocsPositions = false;
    if (bclData.storeXy())
    {
        ISAAC_THREAD_CERR << "Loading Positions data for " << tileMetadata << std::endl;
        flowcell.getLaneAttribute<flowcell::Layout::BclBgzf, flowcell::PositionsFilePathAttributeTag>(tileMetadata.getLane(), positionsFilePath_);
        if (flowcell::isClocsPath(positionsFilePath_))
        {
            ISAAC_ASSERT_MSG(false, "Clocs are not supported with bcl.bgzf data");
            clocsMapper_.mapTile(positionsFilePath_, tileMetadata.getClusterCount());
        }
        else
        {
            locsMapper_.mapTile(positionsFilePath_, tileMetadata.getClusterCount(),
                                tileClusterOffsets_.at(tileBciIndexMap_.at(tileMetadata.getIndex())));
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


void BclBgzfBaseCallsSource::bclToClusters(
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
    ISAAC_ASSERT_MSG(bclData.pf().size() == tileMetadata.getClusterCount(), "Mismatch between expected:" << tileMetadata.getClusterCount() <<
                     " and retrieved:" << bclData.pf().size() << " number of pf values");
    ISAAC_THREAD_CERR << "Extracting Pf values done for " << bclData.getClusterCount() << " bcl clusters" << std::endl;

    if (bclData.storeXy())
    {
        ISAAC_THREAD_CERR << "Extracting Positions values for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
        bclData.xy().clear();
        if (!useLocsPositions)
        {
            ISAAC_ASSERT_MSG(false, "Clocs are not supported with bcl.bgzf data");
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

template class BclBgzfSeedSource<oligo::ShortKmerType>;
template class BclBgzfSeedSource<oligo::KmerType>;
template class BclBgzfSeedSource<oligo::LongKmerType>;

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
