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
 ** \file CycleBciMapper.hh
 **
 ** Helper class for cycle .bci files into memory.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_RTA_CYCLE_BCI_MAPPER_HH
#define iSAAC_RTA_CYCLE_BCI_MAPPER_HH

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"

namespace isaac
{
namespace rta
{

class CycleBciMapper
{
    struct Header
    {
        uint32_t version_;
        uint32_t tilesCount_;
    } __attribute__ ((packed));

public:
    struct VirtualOffset
    {
        uint64_t uncompressedOffset : 16;
        uint64_t compressedOffset : 48;
    } __attribute__ ((packed));

    CycleBciMapper(const std::size_t tilesMax)
    {
        tileOffsets_.reserve(tilesMax);
    }

    CycleBciMapper(const CycleBciMapper &that)
    {
        tileOffsets_.reserve(that.tileOffsets_.capacity());
    }

    void mapFile(const boost::filesystem::path &cycleBciPath)
    {
        std::ifstream is(cycleBciPath.c_str());
        if (!is)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open file: " + cycleBciPath.string() + strerror(errno)));
        }
        mapStream(is, cycleBciPath);
    }

    void mapStream(std::istream &is, const boost::filesystem::path &cycleBciPath)
    {
        Header header;
        if (!is.read(reinterpret_cast<char *>(&header), sizeof(header)))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to reader header from %s: %s") % cycleBciPath % strerror(errno)).str()));
        }

        if (0 != header.version_)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Invalid file format version %d in %s. Expected 0.") % unsigned(header.version_) % cycleBciPath).str()));
        }

        tileOffsets_.resize(header.tilesCount_);
        if (!is.read(reinterpret_cast<char *>(&tileOffsets_.front()), tileOffsets_.size() * sizeof(uint64_t)))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to %d records from %s: %s") % tileOffsets_.size() % cycleBciPath % strerror(errno)).str()));
        }
    }

    /**
     * \brief returns virtual offset in the compressed cycle bcl file given the index of the tile
     * \param tileIndex Index of the tile record in cycle bci. Is is the same as the index of the corresponding tile
     *                  record in lane bci.
     */
    VirtualOffset getTileOffset(const unsigned tileIndex) const
    {
        return tileOffsets_.at(tileIndex);
    }

private:
    std::vector<VirtualOffset> tileOffsets_;
};


} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_RTA_CYCLE_BCI_MAPPER_HH
