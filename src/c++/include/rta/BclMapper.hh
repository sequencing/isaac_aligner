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
 ** \file BclMapper.hh
 **
 ** Helper class for mapping BCL files into memory.
 **
 ** \author Roman Petrovski
 **/
 
#ifndef iSAAC_RTA_BCL_MAPPER_HH
#define iSAAC_RTA_BCL_MAPPER_HH

#include <boost/format.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

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
namespace rta
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
    unsigned long getUnpaddedBclSize() const
    {
        return sizeof(uint32_t) + clusterCount_;
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
    BclMapper(const unsigned maxCycles, const unsigned maxClusters):
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

private:
    typedef boost::error_info<struct tag_errmsg, std::string> errmsg_info;

    unsigned clusterCount_;
    unsigned cycleNumbers_;
    std::vector<char> tileData_;
};

/**
 * \brief helper class to combat the bcl file access latency by loading multiple cycles in parallel
 */
template <typename ReaderT>
class ParallelBclMapper : BclMapper
{
    common::ThreadVector &threads_;
    const unsigned maxInputLoaders_;
    std::vector<ReaderT> &threadReaders_;
    std::vector<unsigned> cycleNumbers_;
public:
    using BclMapper::transpose;
    using BclMapper::getCyclesCount;
    ParallelBclMapper(const bool ignoreMissingBcls,
                      const unsigned maxCycles,
                      common::ThreadVector &threads,
                      std::vector<ReaderT> &threadReaders,
                      const unsigned maxInputLoaders,
                      unsigned maxClusters,
                      const std::size_t reservePathLength):
        BclMapper(maxCycles, maxClusters),
        threads_(threads),
        maxInputLoaders_(maxInputLoaders),
        threadReaders_(threadReaders),
        cycleNumbers_(maxCycles)
    {
        std::for_each(threadReaders_.begin(), threadReaders_.end(),
                      boost::bind(&ReaderT::reserveBuffers, _1, reservePathLength, getTileSize(1)));
    }

    void mapTile(const flowcell::Layout &flowcell, const flowcell::TileMetadata &tileMetadata)
    {
        ISAAC_ASSERT_MSG(cycleNumbers_.capacity() >= flowcell.getDataCycles().size() + flowcell.getBarcodeCycles().size(),
                         "Insufficient capacity in cycleNumbers_ need " << flowcell.getDataCycles().size() + flowcell.getBarcodeCycles().size() << " got " << cycleNumbers_.size());
        cycleNumbers_.clear();

        // Add barcode cycles first
        cycleNumbers_ = flowcell.getBarcodeCycles();
        // Add data cycles second
        cycleNumbers_.insert(cycleNumbers_.end(), flowcell.getDataCycles().begin(), flowcell.getDataCycles().end());

        setGeometry(cycleNumbers_.size(), tileMetadata.getClusterCount());

        threads_.execute(boost::bind(
            &ParallelBclMapper::threadLoadBcls, this, _1,
            tileMetadata.getClusterCount(),
            boost::ref(flowcell), boost::ref(tileMetadata),
            cycleNumbers_.begin(),
            cycleNumbers_.end()), std::min<unsigned>(cycleNumbers_.size(), maxInputLoaders_));
    }

private:
    void threadLoadBcls(
        const unsigned threadNumber,
        const unsigned clusterCount,
        const flowcell::Layout &flowcell, const flowcell::TileMetadata &tileMetadata,
        std::vector<unsigned>::const_iterator threadCyclesBegin,
        std::vector<unsigned>::const_iterator threadCyclesEnd)
    {
        // each thread starts at the threadNumber cycle
        for(unsigned thistThreadCycleOffset = threadNumber;
            std::distance(threadCyclesBegin, threadCyclesEnd) > thistThreadCycleOffset;
            // jump to the next cycle to be loaded by this thread
            thistThreadCycleOffset += threads_.size())
        {
            const unsigned cycle = *(threadCyclesBegin + thistThreadCycleOffset);
            const unsigned readClusters = threadReaders_[threadNumber].readTileCycle(
                flowcell, tileMetadata, cycle,
                getCycleBufferStart(thistThreadCycleOffset), getTileSize(1));
            ISAAC_ASSERT_MSG(readClusters == clusterCount, "Expected Bcl number of clusters(" << clusterCount <<
                             ") does not match the one read from file(readClusters:" << readClusters <<
                             "cycle:" << cycle << "): ");

            //ISAAC_THREAD_CERR << "Read " << clusters << " clusters from " << *threadCyclePathsBegin << std::endl;
        }
    }
};

/**
 * \brief helper class to combat the bcl file access latency by loading multiple cycle s in parallel
 */
template <typename ReaderT>
class SingleCycleBclMapper : BclMapper
{
public:
    using BclMapper::get;
    SingleCycleBclMapper(
        const unsigned maxClusters,
        const std::size_t reservePathLength, const bool reserveCompressionBuffer,
        ReaderT &reader):
        BclMapper(1, maxClusters),
        reader_(reader)
    {
        reader_.reserveBuffers(reservePathLength, reserveCompressionBuffer ? getTileSize(1) : 0);
    }

    void mapTileCycle(const flowcell::Layout &flowcellLayout, const flowcell::TileMetadata &tile, const unsigned cycle)
    {
        setGeometry(1, tile.getClusterCount());

        const unsigned readClusters = reader_.readTileCycle(
            flowcellLayout, tile, cycle, getCycleBufferStart(0), getTileSize(1));

        ISAAC_ASSERT_MSG(readClusters == tile.getClusterCount(), "Expected Bcl number of clusters(" << tile.getClusterCount() <<
                         ") does not match the one read from file(" << readClusters <<
                         "): " << tile << " cycle:" << cycle);
    }

private:
    ReaderT &reader_;
};


} // namespace rta
} // namespace isaac

#endif // #ifndef iSAAC_RTA_BCL_MAPPER_HH
