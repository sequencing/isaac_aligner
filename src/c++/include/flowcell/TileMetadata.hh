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
 ** \file TileMetadata.hh
 **
 ** Packaging of the metadata associated to a tile.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_FLOWCELL_TILE_METADATA_HH
#define iSAAC_FLOWCELL_TILE_METADATA_HH

#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

namespace isaac
{
namespace flowcell
{

/**
 ** \brief Read-only interface to the metadata associated to a tile.
 **
 ** The intended usage is for tile management in ordered collections (the index
 ** in the collectionis is associated to each tile metadata instance).
 **
 ** \todo Provide a robust and flexible indexing mechanism for the tiles.
 **/
class TileMetadata
{
    /*
     * \brief enable serialization
     */
    template <class Archive> friend void serialize(Archive &ar, TileMetadata &tm, const unsigned int version);

public:
    TileMetadata();
    TileMetadata(const TileMetadata &that, const unsigned newIndex);
    TileMetadata(
        const std::string &flowcellId,
        const unsigned flowcellIndex,
        const unsigned tile,
        const unsigned int lane,
        const unsigned int clusterCount,
        const unsigned int index);
    const std::string &getFlowcellId() const {return flowcellId_;}
    unsigned getFlowcellIndex() const {return flowcellIndex_;}
    unsigned int getTile() const {return tile_;}
    const std::string &getTileString() const {return tileString_;}
    unsigned int getLane() const {return lane_;}
    const std::string &getLaneString() const {return laneString_;}

    unsigned int getClusterCount() const {return clusterCount_;}
    unsigned int getIndex() const {return index_;}

private:
    std::string flowcellId_;
    unsigned flowcellIndex_;
    unsigned int tile_;
    std::string tileString_;
    unsigned int lane_;
    std::string laneString_;
    unsigned int clusterCount_;
    unsigned int index_;
};

//typedef std::vector<TileMetadata> TileMetadataList;
struct TileMetadataList : public std::vector<flowcell::TileMetadata>
{
    TileMetadataList(){}
//    TileMetadataList(size_t size) : std::vector<flowcell::TileMetadata>(size){}
    TileMetadataList(const TileMetadataList &that) : std::vector<flowcell::TileMetadata>(that){}
};

inline unsigned getMaxTileClusters(const TileMetadataList &tileMetadataList_)
{
    return tileMetadataList_.empty() ?
        0 : std::max_element(tileMetadataList_.begin(), tileMetadataList_.end(),
                             boost::bind(&flowcell::TileMetadata::getClusterCount, _1)<
                             boost::bind(&flowcell::TileMetadata::getClusterCount, _2))->getClusterCount();
}

inline std::ostream &operator<<(std::ostream &os, const TileMetadata &tileMetadata)
{
    return os << "TileMetadata(" 
              << tileMetadata.getFlowcellId() << ", " 
              << tileMetadata.getTile() << ", " 
              << tileMetadata.getLane() << ", " 
              << tileMetadata.getClusterCount() << ", "
              << tileMetadata.getIndex()
              << ")";
}

inline std::ostream &operator<<(std::ostream &os, const TileMetadataList &tileMetadataList)
{
    std::copy(tileMetadataList.begin(), tileMetadataList.end(), std::ostream_iterator<TileMetadata>(os, " "));
    return os;
}

} // namespace flowcell
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_TILE_METADATA_HH

