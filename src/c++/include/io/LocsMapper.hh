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
 ** \file LocsMapper.hh
 **
 ** Helper class for mapping compressed position files into memory.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_LOCS_MAPPER_HH
#define iSAAC_IO_LOCS_MAPPER_HH

#include "common/Debug.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace io
{

class LocsMapper
{
public:
    LocsMapper() :
        fileBufCache_(1, std::ios_base::in | std::ios_base::binary), clusterCount_(0)
    {
    }

    void mapTile(const boost::filesystem::path &clocsFilePath,
                 const unsigned clusterCount,
                 const unsigned long clusterOffset = ONE_TILE_PER_FILE)
    {
        clusterCount_ = clusterCount;
        tileData_.clear();
        load(clocsFilePath, clusterOffset, V1);
    }

    template <typename InsertIteratorT>
    void getPositions(InsertIteratorT it) const
    {
        getPositions(it, clusterCount_);
    }

    void reserveBuffers(const size_t reservePathLength, const unsigned maxClusterCount)
    {
        fileBufCache_.reservePathBuffers(reservePathLength);
        tileData_.reserve(sizeof(V0Header) + sizeof(Xy) * maxClusterCount);
    }

    void unreserve()
    {
        std::vector<char>().swap(tileData_);
        fileBufCache_.unreserve();
    }
private:
    static const unsigned long ONE_TILE_PER_FILE = -1UL;
    // Bizarre format documented here http://ukch-confluence.illumina.com/display/SWD/RTA+clocs+file+format
    struct V0Header
    {
        uint32_t version_;
        float floatVersion_;
        uint32_t clusters_; //unsigned 32bits little endian integer: number of clusters
    }__attribute__ ((packed));

    struct Xy
    {
        float x_;
        float y_;
    }__attribute__ ((packed));

    io::FileBufCache<io::FileBufWithReopen> fileBufCache_;
    unsigned clusterCount_;
    std::vector<char> tileData_;
    enum Version
    {
        V1 = 1,
    };

    template <typename InsertIteratorT>
    void getPositions(InsertIteratorT it, unsigned clusters) const
    {
        const Xy *currentBlockClusters = reinterpret_cast<const Xy *>(reinterpret_cast<const V0Header *>(&tileData_.front()) + 1);

        while (clusters--)
        {
            *it++ = std::pair<int, int>(
                static_cast<int>(round(1000.0 + 10.0 * currentBlockClusters->x_)),
                static_cast<int>(round(1000.0 + 10.0 * currentBlockClusters->y_)));
            ++currentBlockClusters;
        }
    }

    void load(
        const boost::filesystem::path &locsFilePath,
        const unsigned long clusterOffset,
        Version assumedVersion)
    {
        std::istream is(fileBufCache_.get(locsFilePath));
        if (!is)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to open %s: %s") % locsFilePath % strerror(errno)).str()));
        }

        std::size_t dataSize = sizeof(V0Header) + clusterCount_ * sizeof(Xy);//common::getFileSize(locsFilePath.c_str());
        if (dataSize > tileData_.capacity())
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("locs file is bigger than supported maximum %s: %d. Expected max: ") % locsFilePath % dataSize % tileData_.capacity()).str()));
        }

        tileData_.resize(dataSize);
        if (!is.read(&tileData_.front(), sizeof(V0Header)))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read %d header bytes from %s") %
                sizeof(V0Header) % locsFilePath ).str()));
        }

        V0Header &header = reinterpret_cast<V0Header &>(tileData_.front());
        if (header.version_ != uint32_t(assumedVersion))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Unsupported locs file version %s: %d") % locsFilePath % int(header.version_)).str()));
        }

        if (ONE_TILE_PER_FILE != clusterOffset)
        {
            // patch the clusters number as multitile filter files contain the total number of clusters
            header.clusters_ = clusterCount_;
            const unsigned long clusterByteOffset = clusterOffset * sizeof(Xy);
            if (!is.seekg(clusterByteOffset, is.cur))
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to seek %d bytes in %s") %
                    clusterByteOffset % locsFilePath.string()).str()));
            }
        }

        if (header.clusters_ !=  clusterCount_)
        {
            BOOST_THROW_EXCEPTION(common::IoException(
                errno, (boost::format("Unexpected locs file number of clusters %s: %d. Expected: %d") %
                    locsFilePath % int(header.clusters_) % clusterCount_).str()));
        }

        if (!is.read(&tileData_.front() + sizeof(V0Header), tileData_.size() - sizeof(V0Header)))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read %d bytes from %s") %
                (dataSize - sizeof(V0Header)) % locsFilePath ).str()));
        }

        if (tileData_.size() - sizeof(V0Header) != std::size_t(is.gcount()))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read %d bytes from %s. Read %d") %
                (tileData_.size() - sizeof(V0Header)) % locsFilePath % is.gcount() ).str()));
        }

        ISAAC_THREAD_CERR << "Read " << clusterCount_ << " position values from locs file version " <<
            assumedVersion << ": " << locsFilePath << " cluster offset:" << clusterOffset << std::endl;
    }
};


} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_LOCS_MAPPER_HH
