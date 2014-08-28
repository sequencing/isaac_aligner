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
 ** \file LaneBciMapper.hh
 **
 ** Helper class for lane .bci files into memory.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_RTA_LANE_BCI_MAPPER_HH
#define iSAAC_RTA_LANE_BCI_MAPPER_HH

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/SystemCompatibility.hh"

namespace isaac
{
namespace rta
{

class LaneBciMapper
{
    struct BciRecord
    {
        uint32_t tileNumber_;
        uint32_t tileClusters_;
        friend std::ostream & operator <<(std::ostream &os, const BciRecord& record)
        {
            return os << "BciRecord(" << record.tileNumber_ << "," << record.tileClusters_ << ")";
        }
    } __attribute__ ((packed));
    std::vector<BciRecord> bci_;
    io::FileBufWithReopen fileBuf_;
public:
    LaneBciMapper(const unsigned maxTiles):
        fileBuf_(std::ios_base::binary|std::ios_base::in)
    {
        bci_.reserve(maxTiles);
    }

    void mapFile(const boost::filesystem::path &laneBciPath)
    {
        std::size_t fileSize = common::getFileSize(laneBciPath.c_str());
        ISAAC_THREAD_CERR << "Loading: " << laneBciPath.c_str() << " size:" << fileSize << std::endl;
        std::istream is(fileBuf_.reopen(laneBciPath.c_str(), io::FileBufWithReopen::SequentialOnce));
        if (!is)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open file: " + laneBciPath.string() + strerror(errno)));
        }

        bci_.resize(fileSize / sizeof(BciRecord));
        if (!is.read(reinterpret_cast<char *>(&bci_.front()), fileSize))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to %d bytes from %s: %s") % fileSize % laneBciPath % strerror(errno)).str()));
        }
    }

    unsigned getTilesCount() const
    {
        return bci_.size();
    }

    struct TileInfo
    {
        TileInfo() :
            tileIndex_(MISSING_TILE_INDEX), tileClusters_(0){}
        TileInfo(unsigned tileIndex, unsigned tileClusters):
            tileIndex_(tileIndex), tileClusters_(tileClusters){}

        unsigned tileIndex_;
        unsigned tileClusters_;

        static const unsigned MISSING_TILE_INDEX = -1U;
    };

    /**
     * \return information about the tile or an empty TileInfo if the tile is not present in the bci file
     */
    TileInfo getTileInfo(const unsigned tileNumber) const
    {
        const std::vector<BciRecord>::const_iterator it = std::find_if(
            bci_.begin(), bci_.end(), boost::bind(&BciRecord::tileNumber_, _1) == tileNumber);
        if (bci_.end() == it)
        {
            return TileInfo();
        }
        return TileInfo(std::distance(bci_.begin(), it), it->tileClusters_);
    }

    /**
     * \return number of clusters in the tile given the 0-based index of the tile in lane bci file
     */
    unsigned getTileClusterCount(const unsigned tileIndex) const
    {
        return bci_.at(tileIndex).tileClusters_;
    }

};


} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_RTA_LANE_BCI_MAPPER_HH
