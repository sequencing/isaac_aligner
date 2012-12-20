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
 ** \file BclDataSource.cpp
 **
 ** \brief see BclDataSource.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>

#include "BclDataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{


flowcell::TileMetadataList BclDataSource::discoverTiles()
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

void BclDataSource::initBuffers(
    flowcell::TileMetadataList &unprocessedTiles,
    const alignment::SeedMetadataList &seedMetadataList,
    common::ThreadVector &threads)
{
    seedLoader_.reset(new alignment::ParallelSeedLoader(
        ignoreMissingBcls_, threads,
        inputLoadersMax_, barcodeMetadataList_,
        bclFlowcellLayout_.getReadMetadataList(),
        seedMetadataList,
        sortedReferenceXmlList_, unprocessedTiles));
}

void BclDataSource::generateSeeds(
    const flowcell::TileMetadataList &tiles,
    const alignment::matchFinder::TileClusterInfo &tileClusterBarcode,
    std::vector<alignment::Seed> &seeds,
    common::ScoopedMallocBlock  &mallocBlock)
{
    seedLoader_->loadSeeds(tiles, tileClusterBarcode, seeds, mallocBlock);
}

const std::vector<std::vector<alignment::Seed>::iterator> &BclDataSource::getReferenceSeedBounds()
{
    return seedLoader_->getReferenceSeedBounds();
}

flowcell::TileMetadataList BclDataSource::getTiles(const flowcell::Layout &flowcellLayout) const
{
    flowcell::TileMetadataList tileMetadataList;

    const std::string &flowcellId = flowcellLayout.getFlowcellId();
    const boost::filesystem::path &baseCallsDirectory = flowcellLayout.getBaseCallsDirectory();
    BOOST_FOREACH(const unsigned int lane, flowcellLayout.getLaneIds())
    {
        const std::vector<unsigned int> tileList = flowcellLayout.getTileIds(lane);
        BOOST_FOREACH(const unsigned int tile, tileList)
        {
            const flowcell::TileMetadata::Compression compression = flowcell::Layout::BclGz == flowcellLayout.getFormat()
                    ? flowcell::TileMetadata::GzCompression : flowcell::TileMetadata::NoCompression;
            boost::filesystem::path bclFilePath;
            flowcell::Layout::getBclFilePath(
                tile, lane, baseCallsDirectory,
                flowcellLayout.getAllCycleNumbers().at(0), compression, bclFilePath);
            const unsigned int clusterCount = io::BclMapper::getClusterCount(bclFilePath);
            const flowcell::TileMetadata tileMetadata(
                flowcellId, flowcellLayout.getIndex(),
                tile, lane, baseCallsDirectory,
                clusterCount,
                compression,
                tileMetadataList.size());
            tileMetadataList.push_back(tileMetadata);
        }
    }

    return tileMetadataList;
}


} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
