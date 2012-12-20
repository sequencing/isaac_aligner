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
 ** \file BclMapper.hh
 **
 ** Helper class for mapping BCL files into memory.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_BCL_MAPPER_HH
#define iSAAC_IO_BCL_MAPPER_HH

#include <boost/format.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/pool/pool_alloc.hpp>

#include "common/Debug.hh"
#include "common/FileSystem.hh"
#include "common/Memory.hh"
#include "common/Threads.hpp"
#include "flowcell/Layout.hh"
#include "flowcell/TileMetadata.hh"

#include "io/InflateGzipDecompressor.hh"
#include "io/FileBufCache.hh"

namespace isaac
{
namespace io
{

namespace BclMapperDetails
{
    inline std::ostream &operator << (std::ostream &os, const std::vector<unsigned> &cycleNumbers)
    {
        std::copy(cycleNumbers.begin(), cycleNumbers.end(), std::ostream_iterator<unsigned>(os, ","));
        return os;
    }
}

class BclMapper
{
public:
    template <typename InsertIteratorT>
    void get(unsigned clusterIndex, InsertIteratorT insertIterator) const
    {
        ISAAC_ASSERT_MSG(clusterIndex < clusterCount_, "Requested cluster number is not in the data");
        const char *clusterCycleEnd = getBclBufferStart(cycleNumbers_) + getClusterOffset(clusterIndex);
        const char *clusterCyclePointer = getBclBufferStart(0) + getClusterOffset(clusterIndex);
        const unsigned increment = getTileSize(1);
        while(clusterCycleEnd != clusterCyclePointer)
        {
            *insertIterator++ = *clusterCyclePointer;
            clusterCyclePointer += increment;
        }
    }

    template <typename InsertIteratorT>
    void transpose(InsertIteratorT insertIterator) const
    {
        const unsigned increment = getTileSize(1);
        for (unsigned clusterIndex = 0; clusterCount_ > clusterIndex; ++ clusterIndex)
        {
            const char *clusterCycleEnd = getBclBufferStart(cycleNumbers_) + getClusterOffset(clusterIndex);
            const char *clusterCyclePointer = getBclBufferStart(0) + getClusterOffset(clusterIndex);
            while(clusterCycleEnd != clusterCyclePointer)
            {
                *insertIterator++ = *clusterCyclePointer;
                clusterCyclePointer += increment;
            }
        }
    }

    static unsigned int getClusterCount(const boost::filesystem::path &bclFilePath)
    {
        if (common::isDotGzPath(bclFilePath))
        {
            boost::iostreams::filtering_istream is;
            boost::iostreams::file_source source(bclFilePath.string());
            if (!source.is_open()) {
                BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open file " + bclFilePath.string()));
            }

            is.push(boost::iostreams::gzip_decompressor());
            is.push(source);


            return getClusterCount(is, bclFilePath);
        }
        else
        {
            std::ifstream is(bclFilePath.c_str());
            return getClusterCount(is, bclFilePath);
        }
    }

    void unreserve()
    {
        std::vector<char>().swap(tileData_);
    }

    unsigned getCyclesCount() const {return cycleNumbers_;}

protected:
    void fillWithBlanks(const unsigned cycleIndex)
    {
        char *bufferStart = getCycleBufferStart(cycleIndex);
        std::fill(bufferStart, bufferStart + getTileSize(1), 0);
    }

    unsigned loadFlatBcl(std::istream &is, const unsigned cycleIndex)
    {
        loadRawToTheEnd(is, getCycleBufferStart(cycleIndex), getTileSize(1));
        const unsigned &clusterCount = reinterpret_cast<const unsigned&>(*getCycleBufferStart(cycleIndex));
        ISAAC_ASSERT_MSG(clusterCount == clusterCount_, "Actual Bcl cluster number does not match the one supplied");
        return clusterCount;
    }

    template <typename DecompressorT>
    unsigned loadCompressedBcl(std::istream &source,
                           const unsigned cycleIndex,
                           DecompressorT &decompressor)
    {
        decompressor.reset();
        const unsigned decompressedBytes =
            decompressor.read(source, getCycleBufferStart(cycleIndex), getTileSize(1));

        ISAAC_ASSERT_MSG(decompressedBytes == getUnpaddedBclSize(), "Actual Bcl bytes number does not match the one needed for all clusters");

        const unsigned &clusterCount = reinterpret_cast<const unsigned&>(*getCycleBufferStart(cycleIndex));
        ISAAC_ASSERT_MSG(clusterCount == clusterCount_, "Actual Bcl cluster number does not match the one supplied");
        return clusterCount;
    }

    unsigned loadRawToTheEnd(std::istream &source, char *bufferPos, const unsigned bufferSize)
    {
        using boost::format;
        source.read(bufferPos, bufferSize);
        const unsigned readBytes = source.gcount();

//        source.seekg(0, std::ios::end);
//        ISAAC_ASSERT_MSG(source.tellg() == readBytes, "Expected to reach the end of Bcl stream");
        return readBytes;
    }

    unsigned long getUnpaddedBclSize() const
    {
        return sizeof(clusterCount_) + clusterCount_;
    }

    unsigned long getBclSize() const
    {
        return common::pageRoundUp(getUnpaddedBclSize());
    }

    unsigned long getTileSize(const unsigned cycles) const
    {
        return cycles * getBclSize();
    }

    char *getCycleBufferStart(const unsigned cycleIndex)
    {
        return &tileData_.front() + getTileSize(cycleIndex);
    }

    const char *getBclBufferStart(const unsigned cycleIndex) const
    {
        return &tileData_.front() + getTileSize(cycleIndex);
    }

    const char *getBclBufferEnd(const unsigned cycleIndex) const
    {
        return &tileData_.front() + getTileSize(cycleIndex + 1);
    }

    unsigned getClusterOffset(const unsigned clusterNumber) const
    {
        return sizeof(clusterCount_) + clusterNumber;
    }

    void setGeometry(const unsigned cycles,const unsigned clusterCount)
    {
        clusterCount_ = clusterCount;
        cycleNumbers_ = cycles;
        tileData_.resize(getTileSize(cycleNumbers_));
    }

    /**
     * \brief constructor for child classes that want to perform loading themselves
     */
    BclMapper(const bool ignoreMissingBcls, const unsigned maxCycles, const unsigned maxClusters):
        ignoreMissingBcls_(ignoreMissingBcls),
        clusterCount_(maxClusters),
        cycleNumbers_(maxCycles)
    {
        tileData_.reserve(getTileSize(cycleNumbers_));
    }

    static unsigned int getClusterCount(std::istream &is, const boost::filesystem::path &bclFilePath)
    {
        unsigned int clusterCount;
        char *buffer = reinterpret_cast<char *>(&clusterCount);
        if (!is.read(buffer, 4))
        {
            using boost::format;
            BOOST_THROW_EXCEPTION(std::ios::failure((format("Failed to read cluster count from %s: %s") % bclFilePath % strerror(errno)).str()));
        }
        return clusterCount;
    }

protected:
    const bool ignoreMissingBcls_;
private:
    unsigned clusterCount_;
    unsigned cycleNumbers_;
    std::vector<char> tileData_;
};

/**
 * \brief helper class to combat the bcl file access latency by loading multiple cycle s in parallel
 */
class ParallelBclMapper : BclMapper
{
    common::ThreadVector &threads_;
    boost::ptr_vector<io::InflateGzipDecompressor<std::vector<char> > > decompressors_;
    //std::vector<std::vector<char> > compressedBuffers_;

    std::vector<FileBufCache<FileBufWithReopen> > threadBclFileBuffers_;

public:
    using BclMapper::transpose;
    using BclMapper::getCyclesCount;
    ParallelBclMapper(const bool ignoreMissingBcls,
                      const unsigned maxCycles,
                      common::ThreadVector &threads,
                      unsigned maxClusters):
        BclMapper(ignoreMissingBcls, maxCycles, maxClusters),
        threads_(threads),
        threadBclFileBuffers_(threads_.size(), FileBufCache<FileBufWithReopen>(1, std::ios_base::in | std::ios_base::binary))
    {
    }

    void mapTile(const std::vector<boost::filesystem::path> &cyclePaths, const unsigned clusterCount)
    {
        setGeometry(cyclePaths.size(), clusterCount);

        threads_.execute(boost::bind(
            &ParallelBclMapper::threadLoadBcls, this, _1,
            cyclePaths.begin(),
            cyclePaths.end()));
    }

    void unreserve()
    {
        std::vector<FileBufCache<FileBufWithReopen> >().swap(threadBclFileBuffers_);
        decompressors_.clear();
        BclMapper::unreserve();
    }

    void reserveBuffers(const size_t reservePathLength, const bool reserveCompressionBuffers)
    {
        std::for_each(threadBclFileBuffers_.begin(), threadBclFileBuffers_.end(),
                      boost::bind(&io::FileBufCache<FileBufWithReopen>::reservePathBuffers, _1, reservePathLength));

        if (reserveCompressionBuffers)
        {
            decompressors_.reserve(threads_.size());
            while(decompressors_.size() < threads_.size())
            {
                decompressors_.push_back(new io::InflateGzipDecompressor<std::vector<char> >(getTileSize(1)));
            }
//            compressedBuffers_.resize(threads_.size(), std::vector<char>(getTileSize(1)));
        }
    }

private:
    void threadLoadBcls(
        const unsigned threadNumber,
        std::vector<boost::filesystem::path>::const_iterator threadCyclePathsBegin,
        std::vector<boost::filesystem::path>::const_iterator threadCyclePathsEnd)
    {
        // each thread starts at the threadNumber cycle
        unsigned thistThreadCycle = threadNumber;
        threadCyclePathsBegin += threadNumber;
        while(threadCyclePathsEnd > threadCyclePathsBegin)
        {
            if (ignoreMissingBcls_ && !boost::filesystem::exists(*threadCyclePathsBegin))
            {
                ISAAC_THREAD_CERR << "WARNING: Ignoring missing bcl file: " << *threadCyclePathsBegin << std::endl;
                fillWithBlanks(thistThreadCycle);
            }
            else
            {
                std::istream source(threadBclFileBuffers_[threadNumber].get(*threadCyclePathsBegin, FileBufWithReopen::SequentialOnce));

                /*const unsigned clusters = */common::isDotGzPath(*threadCyclePathsBegin)
                    ? loadCompressedBcl(source, thistThreadCycle, decompressors_.at(threadNumber))
                    : loadFlatBcl(source, thistThreadCycle);
            }

            //ISAAC_THREAD_CERR << "Read " << clusters << " clusters from " << *threadCyclePathsBegin << std::endl;
            // jump to the next cycle to be loaded by this thread
            thistThreadCycle += threads_.size();
            // -D_GLIBCXX_DEBUG=1 does not allow advancing past the  end.
            if (threadCyclePathsEnd - threadCyclePathsBegin <= unsigned(threads_.size()))
            {
                break;
            }
            threadCyclePathsBegin += threads_.size();
        }
    }
};

/**
 * \brief helper class to combat the bcl file access latency by loading multiple cycle s in parallel
 */
class SingleCycleBclMapper : BclMapper
{
public:
    using BclMapper::get;
    SingleCycleBclMapper(const bool ignoreMissingBcls, const unsigned maxClusters):
        BclMapper(ignoreMissingBcls, 1, maxClusters),
        decompressor_(maxClusters),
        bclFileBuffer_(1, std::ios_base::in | std::ios_base::binary)
    {
    }

    void mapTileCycle(const flowcell::TileMetadata &tile, const unsigned cycle)
    {
        setGeometry(1, tile.getClusterCount());

        flowcell::Layout::getBclFilePath(
            tile.getTile(), tile.getLane(), tile.getBaseCallsPath(),
            cycle, tile.getCompression(),
            cycleFilePath_);

        if (ignoreMissingBcls_ && !boost::filesystem::exists(cycleFilePath_))
        {
            ISAAC_THREAD_CERR << "WARNING: Ignoring missing bcl file: " << cycleFilePath_ << std::endl;
            fillWithBlanks(0);
        }
        else
        {
            std::istream source(bclFileBuffer_.get(cycleFilePath_, FileBufWithReopen::SequentialOnce));
            const unsigned clusters = common::isDotGzPath(cycleFilePath_)
                ? loadCompressedBcl(source, 0, decompressor_) : loadFlatBcl(source, 0);
            ISAAC_THREAD_CERR << "Read " << clusters << " clusters from " << cycleFilePath_ << std::endl;
        }
    }

    void reserveBuffers(const size_t reservePathLength, const bool reserveCompressionBuffer)
    {

        bclFileBuffer_.reservePathBuffers(reservePathLength);
        // ensure the cycleFilePath_ owns a buffer of maxFilePathLen capacity
        {cycleFilePath_ = std::string(reservePathLength, 'a');}
        cycleFilePath_.clear();
        if (reserveCompressionBuffer)
        {
            decompressor_.resize(getTileSize(1));
        }
    }

private:
    io::InflateGzipDecompressor<std::vector<char> > decompressor_;
    std::vector<char> compressedBuffer_;

    FileBufCache<FileBufWithReopen> bclFileBuffer_;
    boost::filesystem::path cycleFilePath_;
};


} // namespace io
} // namespace isaac

#endif // #ifndef iSSAC_IO_BCL_MAPPER_HH
