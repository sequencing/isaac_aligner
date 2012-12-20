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
 ** \file BufferingFragmentStorage.cpp
 **
 ** \author Roman Petrovski
 **/

#include <cerrno>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "alignment/BinMetadata.hh"
#include "alignment/matchSelector/BufferingFragmentStorage.hh"
#include "io/Fragment.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

BufferingFragmentStorage::BufferingFragmentStorage(
    const bool keepUnaligned,
    const unsigned maxSavers,
    const unsigned threadBuffers,
    const MatchDistribution &matchDistribution,
    const unsigned long outputBinSize,
    const bfs::path &binDirectory,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const unsigned long maxTileClusters,
    const unsigned long totalTiles)
    : keepUnaligned_(keepUnaligned)
    , maxTileReads_(maxTileClusters * 2)
    , storedTile_(0)
    , binIndexMap_(matchDistribution, outputBinSize)
    , flushThreads_(maxSavers)
    , binPathList_(buildBinPathList(binIndexMap_, matchDistribution.getBinSize(), binDirectory, barcodeMetadataList, maxTileReads_, totalTiles))
    , fragmentCollector_(binIndexMap_, maxTileClusters, flowcellLayoutList)
    , flushBuffer_(maxTileClusters, flowcellLayoutList)
    , threadDataFileBufCaches_(flushThreads_.size(),
                           FileBufCache(1, std::ios_base::out | std::ios_base::app | std::ios_base::binary))
    , threadFIdxFileBufCaches_(flushThreads_.size(),
                       FileBufCache(1, std::ios_base::out | std::ios_base::app | std::ios_base::binary))
    , threadRIdxFileBufCaches_(flushThreads_.size(),
                       FileBufCache(1, std::ios_base::out | std::ios_base::app | std::ios_base::binary))
    , threadSeIdxFileBufCaches_(flushThreads_.size(),
                       FileBufCache(1, std::ios_base::out | std::ios_base::app | std::ios_base::binary))

{
    ISAAC_THREAD_CERR << "Resetting ouput files for " << binIndexMap_.getBinCount() << " bins" << std::endl;

    BOOST_FOREACH(const BinMetadata &binMetadata, binPathList_)
    {
        // Create or truncate the destination file
        std::ofstream os(binMetadata.getPathString().c_str());
        if (!os)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open bin file " + binMetadata.getPathString()));
        }

        // Create or truncate the destination f-idx file
        std::ofstream osFIdx(binMetadata.getFIdxFilePath().c_str());
        if (!osFIdx)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open bin index file " + binMetadata.getFIdxFilePath().string()));
        }

        // Create or truncate the destination r-idx file
        std::ofstream osRIdx(binMetadata.getRIdxFilePath().c_str());
        if (!osRIdx)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open bin index file " + binMetadata.getRIdxFilePath().string()));
        }

        // Create or truncate the destination se-idx file
        std::ofstream osSeIdx(binMetadata.getSeIdxFilePath().c_str());
        if (!osSeIdx)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open bin index file " + binMetadata.getSeIdxFilePath().string()));
        }
    }

    // assuming the last entry in the list contains the longest paths
    std::for_each(threadDataFileBufCaches_.begin(), threadDataFileBufCaches_.end(),
                  boost::bind(&FileBufCache::reservePathBuffers, _1,
                              binPathList_.back().getPathString().size()));
    std::for_each(threadFIdxFileBufCaches_.begin(), threadFIdxFileBufCaches_.end(),
                  boost::bind(&FileBufCache::reservePathBuffers, _1,
                              binPathList_.back().getFIdxFilePath().string().size()));
    std::for_each(threadRIdxFileBufCaches_.begin(), threadRIdxFileBufCaches_.end(),
                  boost::bind(&FileBufCache::reservePathBuffers, _1,
                              binPathList_.back().getRIdxFilePath().string().size()));
    std::for_each(threadSeIdxFileBufCaches_.begin(), threadSeIdxFileBufCaches_.end(),
                  boost::bind(&FileBufCache::reservePathBuffers, _1,
                              binPathList_.back().getSeIdxFilePath().string().size()));


//    threadFragmentCollectors_.reserve(reserveClusters);

    ISAAC_THREAD_CERR << "Resetting ouput files done for " << binIndexMap_.getBinCount() << " bins" << std::endl;
}

void BufferingFragmentStorage::flushBin(
    std::ostream &osData, std::ostream &osFIdx, std::ostream &osRIdx, std::ostream &osSeIdx,
    const unsigned binNumber, BinMetadata &binMetadata)
{
//    ISAAC_THREAD_CERR_DEV_TRACE((boost::format("BufferingFragmentStorage::flushBin: %d %s") % binNumber % binMetadata).str());

    ISAAC_ASSERT_MSG(binNumber == binMetadata.getIndex(), "Bin index mismatch");

//    ISAAC_THREAD_CERR_DEV_TRACE((boost::format("flushBuffer_.indexEnd() - currentBinIterator: %d") %
//        (flushBuffer_.indexEnd() - currentBinIterator)).str());

    const FragmentBuffer::IndexIterator binBegin = flushBuffer_.binBegin(binNumber, binIndexMap_);
    FragmentBuffer::IndexIterator binEnd = flushBuffer_.binBegin(binNumber + 1, binIndexMap_);
    FragmentBuffer::IndexIterator currentBinIterator = binBegin;
    // update the offsets and mate offsets assuming all data will be loaded sequentially from the file
    // this has to be a separate pass as reads further down will update their preceding mates
    for(;
        binEnd != currentBinIterator && currentBinIterator->initialized();
        ++currentBinIterator)
    {
        FragmentBuffer::IndexRecord &recordStart = *currentBinIterator;
        const io::FragmentHeader& header = recordStart.fragmentHeader();

        currentBinIterator->setDataOffset(
            binMetadata.incrementDataSize(
                header.fStrandPosition_,
                header.getTotalLength()).first,
            header.flags_.paired_ &&
            binNumber == recordStart.peIndex().mate_.info_.fields_.storageBin_);
    }
    binEnd = currentBinIterator;

    // store data sequentially in the bin file
    for(FragmentBuffer::IndexConstIterator currentBinIterator = binBegin;
        binEnd != currentBinIterator;
        ++currentBinIterator)
    {
        const FragmentBuffer::IndexRecord &recordStart = *currentBinIterator;

//        ISAAC_THREAD_CERR_DEV_TRACE((boost::format("%s") % recordStart).str());

        const io::FragmentHeader& header = recordStart.fragmentHeader();
//        ISAAC_ASSERT_MSG(io::FragmentHeader::magicValue_ == header.magic_, "corrupt binary data in memory");
//        ISAAC_ASSERT_MSG(header.getTotalLength() == header.totalLength_, "corrupt binary data in memory. fragment total length is bad");

        if (!osData.write(reinterpret_cast<const char *>(&header), header.getTotalLength())) {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write into " + binMetadata.getPathString()));
        }

        if (!header.flags_.paired_)
        {
            const io::SeFragmentIndex &idx = recordStart.seIndex();

            if (!osSeIdx.write(reinterpret_cast<const char*>(&idx), sizeof(idx))) {
                BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write into " + binMetadata.getSeIdxFilePath().string()));
            }
            binMetadata.incrementSeIdxElements(header.fStrandPosition_, 1, header.barcode_);
        }
        else if (header.flags_.reverse_ || header.flags_.unmapped_)
        {
//                        ISAAC_THREAD_CERR << "r cluster " << clusterId << " bin" << binNumber << std::endl;
            const io::RStrandOrShadowFragmentIndex &idx = recordStart.rsIndex();

            if (!osRIdx.write(reinterpret_cast<const char*>(&idx), sizeof(idx))) {
                BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write into " + binMetadata.getRIdxFilePath().string()));
            }
            binMetadata.incrementRIdxElements(header.fStrandPosition_, 1, header.barcode_);
        }
        else
        {
//                        ISAAC_THREAD_CERR << "f cluster " << clusterId << " bin" << binNumber << std::endl;
            const io::FStrandFragmentIndex &idx = recordStart.fIndex();

            if (!osFIdx.write(reinterpret_cast<const char*>(&idx), sizeof(idx))) {
                BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write into " + binMetadata.getFIdxFilePath().string()));
            }
            binMetadata.incrementFIdxElements(header.fStrandPosition_, 1, header.barcode_);
        }
        binMetadata.incrementGapCount(header.fStrandPosition_, header.gapCount_, header.barcode_);
        binMetadata.incrementCigarLength(header.fStrandPosition_, header.cigarLength_, header.barcode_);
    }
//    ISAAC_THREAD_CERR << "Flushing bin done: " << binMetadata << std::endl;
}


void BufferingFragmentStorage::flushUnmappedBin(
    std::ostream &osData,
    const unsigned binNumber, BinMetadata &binMetadata)
{
//    ISAAC_THREAD_CERR_DEV_TRACE((boost::format("BufferingFragmentStorage::flushUnmappedBin: %d %s") % binNumber % binMetadata).str());

    ISAAC_ASSERT_MSG(binNumber == binMetadata.getIndex(), "Bin index mismatch");
    ISAAC_ASSERT_MSG(0 == binNumber, "Unaligned in index must be 0");

    const FragmentBuffer::IndexIterator binBegin = flushBuffer_.binBegin(binNumber, binIndexMap_);
    FragmentBuffer::IndexIterator binEnd = flushBuffer_.binBegin(binNumber + 1, binIndexMap_);
    FragmentBuffer::IndexIterator currentBinIterator = binBegin;
    unsigned long storedTileRead=0;
    for(;
        binEnd != currentBinIterator && currentBinIterator->initialized();
        ++currentBinIterator)
    {
        FragmentBuffer::IndexRecord &recordStart = *currentBinIterator;
        const io::FragmentAccessor& fragment = recordStart.fragment();

        currentBinIterator->setDataOffset(binMetadata.getDataSize(), fragment.flags_.paired_);
        binMetadata.incrementDataSize(maxTileReads_ * storedTile_ + storedTileRead++, fragment.getTotalLength());
    }

    binEnd = currentBinIterator;
    storedTileRead = 0;
    for(FragmentBuffer::IndexConstIterator currentBinIterator = binBegin;
        binEnd != currentBinIterator;
        ++currentBinIterator)
    {
        const FragmentBuffer::IndexRecord &recordStart = *currentBinIterator;

        const io::FragmentAccessor& fragment = recordStart.fragment();

        if (!osData.write(reinterpret_cast<const char *>(&fragment), fragment.getTotalLength())) {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write into " + binMetadata.getPathString()));
        }
        binMetadata.incrementNmElements(maxTileReads_ * storedTile_ + storedTileRead++, 1, fragment.barcode_);
    }
//    ISAAC_THREAD_CERR << "Flushing bin done: " << binMetadata << std::endl;
}


void BufferingFragmentStorage::threadFlushBins(const unsigned threadNumber, unsigned &nextUnflushedBin)
{
    boost::lock_guard<boost::mutex> lock(binFlushMutex_);
    while(binIndexMap_.getBinCount() > nextUnflushedBin)
    {
        const unsigned ourBin = nextUnflushedBin++;

        {
            common::unlock_guard<boost::mutex > unlock(binFlushMutex_);

            std::ostream osData(threadDataFileBufCaches_[threadNumber].get(binPathList_[ourBin].getPath()));
            if (0 == ourBin)
            {
                if (keepUnaligned_)
                {
                    flushUnmappedBin(osData, ourBin, binPathList_[ourBin]);
                }
            }
            else
            {
                std::ostream osFIdx(threadFIdxFileBufCaches_[threadNumber].get(binPathList_[ourBin].getFIdxFilePath()));
                std::ostream osRIdx(threadRIdxFileBufCaches_[threadNumber].get(binPathList_[ourBin].getRIdxFilePath()));
                std::ostream osSeIdx(threadSeIdxFileBufCaches_[threadNumber].get(binPathList_[ourBin].getSeIdxFilePath()));
                flushBin(osData, osFIdx, osRIdx, osSeIdx, ourBin, binPathList_[ourBin]);
                // it is very important to flush these files after we're done. Although the FileBufWithReopen will do pubsync
                // itself, most likely the same file will be reopen on a different thread. this means, another FileBufWithReopen
                // will be already writing to it.
                osFIdx.flush();
                osRIdx.flush();
                osSeIdx.flush();
            }
            // it is very important to flush these files after we're done. Although the FileBufWithReopen will do pubsync
            // itself, most likely the same file will be reopen on a different thread. this means, another FileBufWithReopen
            // will be already writing to it.
            osData.flush();
        }
    }
}

void BufferingFragmentStorage::prepareFlush()
{
    ISAAC_ASSERT_MSG(flushBuffer_.empty(), "Buffer must be flushed at this point");
    fragmentCollector_.swapBuffer(flushBuffer_);
}

void BufferingFragmentStorage::flush()
{
    ISAAC_THREAD_CERR << "Flushing buffer" << std::endl;

    flushBuffer_.sortIndex(binIndexMap_);

    unsigned nextUnflushedBin = 0;
    flushThreads_.execute(boost::bind(
            &BufferingFragmentStorage::threadFlushBins, this, _1,
            boost::ref(nextUnflushedBin)));

    assert(binIndexMap_.getBinCount() == nextUnflushedBin);

    flushBuffer_.clear();
    ++storedTile_;

    ISAAC_THREAD_CERR << "Flushing buffer done for " << nextUnflushedBin << " bins" << std::endl;
}

alignment::BinMetadataList BufferingFragmentStorage::buildBinPathList(
    const BinIndexMap &binIndexMap,
    const unsigned long outputBinSize,
    const bfs::path &binDirectory,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const unsigned long maxTileReads,
    const unsigned long totalTiles)
{
    ISAAC_THREAD_CERR << "maxTileClusters " << maxTileReads << "totalTiles " << totalTiles << std::endl;
    alignment::BinMetadataList binPathList;
    ISAAC_ASSERT_MSG(!binIndexMap.empty(), "Empty binIndexMap is illegal");
    ISAAC_ASSERT_MSG(!binIndexMap.back().empty(), "Empty binIndexMap entry is illegal");
    binPathList.reserve(1 + binIndexMap.back().back());
    size_t contigIndex = 0;
    BOOST_FOREACH(const std::vector<unsigned> &contigBins, binIndexMap)
    {
        assert(!contigBins.empty());
        for (unsigned i = contigBins.front(); contigBins.back() >= i; ++i)
        {
            ISAAC_ASSERT_MSG(binPathList.size() == i, "Basic sanity checking for bin numbering failed");
            using boost::format;
            // padd file names well, so that we don't have to worry about them becoming of different length.
            // This is important for memory reservation to be stable
            const reference::ReferencePosition binStartPos = binIndexMap.getBinFirstPos(i);
            binPathList.push_back(
                alignment::BinMetadata(
                    barcodeMetadataList.size(),
                    binPathList.size(),
                    binStartPos,
                    // bin zero has length of -1U as it contains unaligned records which are chunked by 32 bases of their sequence
                    i ? binIndexMap.getBinFirstInvalidPos(i) - binStartPos : maxTileReads * totalTiles,
                    binDirectory / (format("bin-%04d-%04d.dat") % contigIndex % i).str()));
        }
        ++contigIndex;
    }
    return binPathList;
}

} //namespace matchSelector
} // namespace alignment
} // namespace isaac
