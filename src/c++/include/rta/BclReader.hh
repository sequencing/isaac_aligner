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
 ** \file BclReader.hh
 **
 ** Helper class for reading flat and compressed BCL files into memory.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_RTA_BCL_GZ_READER_HH
#define iSAAC_RTA_BCL_GZ_READER_HH

#include <boost/format.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "common/Debug.hh"
#include "common/FileSystem.hh"
#include "common/Memory.hh"
#include "common/Threads.hpp"
#include "flowcell/BclLayout.hh"
#include "flowcell/TileMetadata.hh"

#include "io/InflateGzipDecompressor.hh"
#include "io/FileBufCache.hh"

namespace isaac
{
namespace rta
{

class BclReader
{
public:
    BclReader(const BclReader &that) :
        ignoreMissingBcls_(that.ignoreMissingBcls_),
        decompressor_(that.decompressor_),
        bclFileBuffer_(1, std::ios_base::in | std::ios_base::binary)
    {
    }

    BclReader(const bool ignoreMissingBcls, const unsigned maxClusters):
        ignoreMissingBcls_(ignoreMissingBcls),
        decompressor_(maxClusters),
        bclFileBuffer_(1, std::ios_base::in | std::ios_base::binary)
    {
    }

    void reserveBuffers(
        const std::size_t reservePathLength,
        const std::size_t maxDecompressedBytes)
    {
        // ensure the cycleFilePath_ owns a buffer of maxFilePathLen capacity
        {cycleFilePath_ = std::string(reservePathLength, 'a');}
        cycleFilePath_.clear();

        bclFileBuffer_.reservePathBuffers(reservePathLength);
        if (maxDecompressedBytes)
        {
            decompressor_.resize(maxDecompressedBytes);
        }
    }

    unsigned readTileCycle(
        const flowcell::Layout &flowcellLayout,
        const flowcell::TileMetadata &tile,
        const unsigned cycle,
        char *cycleBuffer, const std::size_t bufferSize)
    {
        flowcellLayout.getLaneTileCycleAttribute<flowcell::Layout::Bcl, flowcell::BclFilePathAttributeTag>(
            tile.getLane(), tile.getTile(), cycle, cycleFilePath_);

        if (ignoreMissingBcls_ && !boost::filesystem::exists(cycleFilePath_))
        {
            ISAAC_THREAD_CERR << "WARNING: Ignoring missing bcl file: " << cycleFilePath_ << std::endl;
            std::fill(cycleBuffer, cycleBuffer + tile.getClusterCount(), 0);
            return tile.getClusterCount();
        }
        else
        {
            std::istream source(bclFileBuffer_.get(cycleFilePath_, io::FileBufWithReopen::SequentialOnce));
            const unsigned clusters = common::isDotGzPath(cycleFilePath_) ?
                loadCompressedBcl(source, cycleFilePath_, cycleBuffer, bufferSize) :
                loadFlatBcl(source, cycleFilePath_, cycleBuffer, bufferSize);
//            ISAAC_THREAD_CERR << "Read " << clusters << " clusters from " << cycleFilePath_ << std::endl;
            return clusters;
        }
    }

private:
    const bool ignoreMissingBcls_;
    boost::filesystem::path cycleFilePath_;
    io::InflateGzipDecompressor<std::vector<char> > decompressor_;

    io::FileBufCache<io::FileBufWithReopen> bclFileBuffer_;

    typedef boost::error_info<struct tag_errmsg, std::string> errmsg_info;

    unsigned loadFlatBcl(std::istream &is,
                         const boost::filesystem::path &filePath,
                         char *bufferStart, const std::size_t bufferSize)
    {
        loadRawToTheEnd(is, bufferStart, bufferSize);
        const unsigned &clusterCount = reinterpret_cast<const unsigned&>(*bufferStart);
        return clusterCount;
    }

    unsigned loadCompressedBcl(std::istream &source,
                           const boost::filesystem::path &filePath,
                           char *bufferStart,
                           const std::size_t bufferSize)
    {
        try
        {
            decompressor_.reset();
            const unsigned decompressedBytes = decompressor_.read(source, bufferStart, bufferSize);

            ISAAC_ASSERT_MSG(4 <= decompressedBytes,
                             "Size of uncompressed bcl data is less than absolute minimum. Required >=4, got:" <<
                             decompressedBytes)
            const unsigned &clusterCount = reinterpret_cast<const unsigned&>(*bufferStart);
            ISAAC_ASSERT_MSG(sizeof(clusterCount) + clusterCount == decompressedBytes,
                             "Actual Bcl bytes number (" << decompressedBytes <<
                             ") does not match the one needed for all clusters (" << (sizeof(clusterCount) + clusterCount) <<
                             ") in file: " << filePath);

            return clusterCount;
        }
        catch (boost::exception &e)
        {
            e << errmsg_info(" While reading from " + filePath.string());
            throw;
        }
        return 0;
    }

    unsigned loadRawToTheEnd(std::istream &source, char *bufferPos, const unsigned bufferSize)
    {
        using boost::format;
        source.read(bufferPos, bufferSize);
        const unsigned readBytes = source.gcount();
        return readBytes;
    }

};


} // namespace rta
} // namespace isaac

#endif // #ifndef iSAAC_RTA_BCL_GZ_READER_HH
