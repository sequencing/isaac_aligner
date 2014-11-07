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
 ** \file BclBgzfTileReader.hh
 **
 ** Helper class for mapping BCL files into memory.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_RTA_BCL_BGZF_TILE_READER_HH
#define iSAAC_RTA_BCL_BGZF_TILE_READER_HH

#include <boost/format.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "common/Debug.hh"
#include "common/FileSystem.hh"
#include "common/Memory.hh"
#include "common/Threads.hpp"
#include "flowcell/BclBgzfLayout.hh"
#include "flowcell/TileMetadata.hh"
#include "io/InflateGzipDecompressor.hh"
#include "io/FileBufCache.hh"
#include "rta/CycleBciMapper.hh"

namespace isaac
{
namespace rta
{

class CycleTileOffsetMap
{

};


class BclBgzfTileReader
{
public:
    BclBgzfTileReader(const BclBgzfTileReader &that) :
        ignoreMissingBcls_(that.ignoreMissingBcls_),
        tileBciIndexMap_(that.tileBciIndexMap_),
        cycleBciMappers_(that.cycleBciMappers_),
        decompressor_(that.decompressor_),
        bclFileBuffer_(std::ios_base::in | std::ios_base::binary)
    {
    }

    BclBgzfTileReader(
        const bool ignoreMissingBcls,
        const unsigned maxClusters,
        const std::vector<unsigned> &tileBciIndexMap,
        const std::vector<rta::CycleBciMapper> &cycleBciMappers):
        ignoreMissingBcls_(ignoreMissingBcls),
        tileBciIndexMap_(tileBciIndexMap),
        cycleBciMappers_(cycleBciMappers),
        decompressor_(maxClusters),
        bclFileBuffer_(std::ios_base::in | std::ios_base::binary)
    {
    }

    void reserveBuffers(
        const std::size_t reservePathLength,
        const std::size_t maxDecompressedBytes)
    {
        // ensure the cycleFilePath_ owns a buffer of maxFilePathLen capacity
        {cycleFilePath_ = std::string(reservePathLength, 'a');}
        cycleFilePath_.clear();

        openFilePath_ = std::string(reservePathLength, 'a');
        decompressor_.resize(maxDecompressedBytes);
    }

    unsigned readTileCycle(
        const flowcell::Layout &flowcellLayout,
        const flowcell::TileMetadata &tile,
        const unsigned cycle,
        char *cycleBuffer, const std::size_t bufferSize)
    {
        flowcellLayout.getLaneCycleAttribute<
            flowcell::Layout::BclBgzf,flowcell::BclFilePathAttributeTag>(tile.getLane(), cycle, cycleFilePath_);
        ISAAC_ASSERT_MSG(tile.getClusterCount() < bufferSize, "Insufficient buffer to read all clusters for " <<
                         cycleFilePath_ << " bufferSize:" << bufferSize << tile);

        if (ignoreMissingBcls_ && !boost::filesystem::exists(cycleFilePath_))
        {
            ISAAC_THREAD_CERR << "WARNING: Ignoring missing bcl file: " << cycleFilePath_ << std::endl;
            std::fill(cycleBuffer, cycleBuffer + tile.getClusterCount(), 0);
            return tile.getClusterCount();
        }
        else
        {
            std::istream source(
                openFilePath_ == cycleFilePath_ ? &bclFileBuffer_ :
                    // keep the cycle file open as we're continuing to read the next tile from the same file
                    bclFileBuffer_.reopen(cycleFilePath_.c_str(), io::FileBufWithReopen::SequentialOnce));
            openFilePath_ = cycleFilePath_.c_str(); // avoid string buffer sharing on copy
            *reinterpret_cast<boost::uint32_t*>(cycleBuffer) = tile.getClusterCount();

            const unsigned clusters = loadCompressedBcl(
                source, cycleFilePath_,
                cycleBciMappers_.at(cycle).getTileOffset(tileBciIndexMap_.at(tile.getIndex())),
                cycleBuffer + sizeof(boost::uint32_t), tile.getClusterCount());
//            ISAAC_THREAD_CERR << "Read " << clusters << " clusters from " << cycleFilePath_ << std::endl;
            return clusters;
        }
    }

private:
    const bool ignoreMissingBcls_;
    const std::vector<unsigned> &tileBciIndexMap_;
    const std::vector<rta::CycleBciMapper> &cycleBciMappers_;
    io::InflateGzipDecompressor<std::vector<char> > decompressor_;
    boost::filesystem::path cycleFilePath_;
    io::FileBufWithReopen bclFileBuffer_;
    boost::filesystem::path openFilePath_;

    typedef boost::error_info<struct tag_errmsg, std::string> errmsg_info;

    unsigned loadCompressedBcl(std::istream &source,
                           const boost::filesystem::path &filePath,
                           const rta::CycleBciMapper::VirtualOffset &tileOffset,
                           char *bufferStart,
                           const std::size_t bufferSize)
    {
        try
        {
            decompressor_.reset();
            if (!source.seekg(tileOffset.compressedOffset))
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to seek to position %d in %s: %s") %
                    boost::uint64_t(tileOffset.compressedOffset) % filePath % strerror(errno)).str()));
            }
            const unsigned decompressedBytes = decompressor_.read(source, tileOffset.uncompressedOffset, bufferStart, bufferSize);

            return decompressedBytes;
        }
        catch (boost::exception &e)
        {
            e << errmsg_info(" While reading from " + filePath.string());
            throw;
        }
        return 0;
    }
};

} // namespace rta
} // namespace isaac

#endif // #ifndef iSAAC_RTA_BCL_BGZF_TILE_READER_HH
