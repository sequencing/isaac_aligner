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
 ** \file TileMetadata.cpp
 **
 ** Packaging of the metadata associated to a tile.
 **
 ** \author Come Raczy
 **/

#include <fstream>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace flowcell
{

TileMetadata::TileMetadata()
    : flowcellId_()
    , flowcellIndex_(-1U)
    , tile_(-1U)
    , tileString_()
    , lane_(-1U)
    , laneString_()
    , clusterCount_(-1U)
    , index_(-1U)
{
}

TileMetadata::TileMetadata(
    const std::string &flowcellId,
    const unsigned flowcellIndex,
    const unsigned tile,
    const unsigned int lane,
    const unsigned int clusterCount,
    const unsigned int index)
    : flowcellId_(flowcellId)
    , flowcellIndex_(flowcellIndex)
    , tile_(tile)
    , tileString_(boost::lexical_cast<std::string>(tile))
    , lane_(lane)
    , laneString_(boost::lexical_cast<std::string>(lane))
    , clusterCount_(clusterCount)
    , index_(index)
{
}

TileMetadata::TileMetadata(const TileMetadata &tileMetadata)
    : flowcellId_(tileMetadata.flowcellId_)
    , flowcellIndex_(tileMetadata.flowcellIndex_)
    , tile_(tileMetadata.tile_)
    , tileString_(boost::lexical_cast<std::string>(tileMetadata.tile_))
    , lane_(tileMetadata.lane_)
    , laneString_(boost::lexical_cast<std::string>(tileMetadata.lane_))
    , clusterCount_(tileMetadata.clusterCount_)
    , index_(tileMetadata.index_)
{
}

TileMetadata::TileMetadata(const TileMetadata &that, const unsigned newIndex)
{
    *this = that;
    index_ = newIndex;
}


TileMetadata &TileMetadata::operator=(const TileMetadata &tileMetadata)
{
    if (this != &tileMetadata)
    {
        setFlowcellId(tileMetadata.flowcellId_);
        flowcellIndex_ = tileMetadata.flowcellIndex_;
        setTile(tileMetadata.tile_);
        setLane(tileMetadata.lane_);
        setClusterCount(tileMetadata.clusterCount_);
        setIndex(tileMetadata.index_);
    }
    return *this;
}

void TileMetadata::setFlowcellId(const std::string &flowcellId)
{
    flowcellId_ = flowcellId;
}

void TileMetadata::setTile(const unsigned int tile)
{
    tile_ = tile;
    tileString_ = boost::lexical_cast<std::string>(tile_);
}

void TileMetadata::setLane(const unsigned int lane)
{
    lane_ = lane;
    laneString_ = boost::lexical_cast<std::string>(lane_);
}

void TileMetadata::setClusterCount(const unsigned int clusterCount)
{
    clusterCount_ = clusterCount;
}

void TileMetadata::setIndex(const unsigned int index)
{
    index_ = index;
}

bool TileMetadata::operator==(const TileMetadata &rhs) const
{
    return
        flowcellId_ == rhs.flowcellId_ &&
        flowcellIndex_ == rhs.flowcellIndex_ &&
        tile_ == rhs.tile_ &&
        lane_ == rhs.lane_ &&
        clusterCount_ == rhs.clusterCount_ &&
        index_ == rhs.index_;
}

} // namespace flowcell
} // namespace isaac
