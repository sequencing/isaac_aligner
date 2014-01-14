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
 ** \file ClocsMapper.hh
 **
 ** Helper class for mapping compressed position files into memory.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_CLOCS_MAPPER_HH
#define iSAAC_IO_CLOCS_MAPPER_HH

#include "common/Debug.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace io
{

class ClocsMapper
{
public:
    ClocsMapper() :
        fileBufCache_(1, std::ios_base::in | std::ios_base::binary), clusterCount_(0)
    {
    }

    void mapTile(const boost::filesystem::path &clocsFilePath, const unsigned clusterCount)
    {
        clusterCount_ = clusterCount;
        tileData_.clear();
        load(clocsFilePath, V1);
    }

    template <typename InsertIteratorT>
    void getPositions(InsertIteratorT it) const
    {
        getPositions(it, clusterCount_);
    }

    void reserveBuffers(const size_t reservePathLength, const unsigned maxClusterCount)
    {
        fileBufCache_.reservePathBuffers(reservePathLength);
        tileData_.reserve(FILE_BYTES_MAX);
    }

    void unreserve()
    {
        std::vector<char>().swap(tileData_);
        fileBufCache_.unreserve();
    }
private:
    // Bizarre format documented here http://ukch-confluence.illumina.com/display/SWD/RTA+clocs+file+format
    struct V0Header
    {
        struct Header
        {
            char version_;
            unsigned bocks_; //unsigned 32bits little endian integer: number of blocks
        }__attribute__ ((packed)) header_;

        struct Block
        {
            unsigned char clusters_;
            struct BlockOffset
            {
                unsigned char blockX_;
                unsigned char blockY_;
            }__attribute__ ((packed)) xy_[0];
        }__attribute__ ((packed)) blocks_[1];
    }__attribute__ ((packed));

    static const int BLOCK_SIZE_X = 25;
    static const int BLOCK_SIZE_Y = 25;
    static const int IMAGE_WIDTH = 2048;
    static const int BLOCKS_PER_LINE = (IMAGE_WIDTH + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X;
    static const int IMAGE_HEIGHT = 20000;
    static const int BLOCKS_PER_COLUMN = (IMAGE_HEIGHT + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y;
    static const int BLOCK_BYTES_MAX = 1 + 255 * 2;// count of clusters plus max number of clusters times two bytes
    static const std::size_t FILE_BYTES_MAX = sizeof(V0Header::Header) + BLOCKS_PER_LINE * BLOCKS_PER_COLUMN * BLOCK_BYTES_MAX;

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
        const V0Header &header = reinterpret_cast<const V0Header &>(tileData_.front());

        const V0Header::Block *currentBlock = header.blocks_;
        int currentBlockX = 0;
        int currentBlockY = 0;
        while (clusters)
        {
            ISAAC_ASSERT_MSG(reinterpret_cast<const char *>(currentBlock) <= &tileData_.back(), "Went outside the clocs file content.");
            unsigned char currentBlockClusters = currentBlock->clusters_;
            const V0Header::Block::BlockOffset *currentBlockCluster = currentBlock->xy_;
            while (currentBlockClusters--)
            {
                *it++ = std::pair<int, int>(currentBlockX + currentBlockCluster->blockX_, currentBlockY + currentBlockCluster->blockY_);
                ++currentBlockCluster;
                ISAAC_ASSERT_MSG(clusters, "More clusters in the clocs file than described in the header");
            }
            currentBlockX += BLOCK_SIZE_X;
            if (IMAGE_WIDTH <= currentBlockX)
            {
                currentBlockX = 0;
                currentBlockY += BLOCK_SIZE_Y;
            }
            const std::size_t lastBlockClusters = std::distance(currentBlock->xy_, currentBlockCluster);
            clusters -= lastBlockClusters;
            currentBlock = reinterpret_cast<const V0Header::Block *>(reinterpret_cast<const char *>(currentBlock) +
                sizeof(V0Header::Block::clusters_) + lastBlockClusters * sizeof(V0Header::Block::BlockOffset));
        }
    }

    void load(const boost::filesystem::path &clocsFilePath, Version assumedVersion)
    {
        std::istream is(fileBufCache_.get(clocsFilePath));
        if (!is)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to open %s: %s") % clocsFilePath % strerror(errno)).str()));
        }

        std::size_t fileSize = common::getFileSize(clocsFilePath.c_str());
        if (fileSize > FILE_BYTES_MAX)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Clocs file is bigger than supported maximum %s: %d") % clocsFilePath % fileSize).str()));
        }

        tileData_.resize(fileSize);
        if (!is.read(&tileData_.front(), fileSize))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read %d bytes from %s") %
                fileSize % clocsFilePath ).str()));
        }

        const V0Header &header = reinterpret_cast<const V0Header &>(tileData_.front());
        if (header.header_.version_ !=  assumedVersion)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Unsupported clocs file version %s: %d") % clocsFilePath % int(header.header_.version_)).str()));
        }

        ISAAC_THREAD_CERR << "Read " << clusterCount_ << " position values from clocs file version " << assumedVersion << ": " << clocsFilePath << std::endl;
    }
};


} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_CLOCS_MAPPER_HH
