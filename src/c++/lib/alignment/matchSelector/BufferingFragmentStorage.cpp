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
    const bool preSortBins,
    const unsigned maxSavers,
    const unsigned threadBuffers,
    const MatchDistribution &matchDistribution,
    const unsigned long outputBinSize,
    const bfs::path &binDirectory,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const unsigned long maxTileClusters,
    const unsigned long totalTiles,
    const bool skipEmptyBins)
    : keepUnaligned_(keepUnaligned)
    , maxTileReads_(maxTileClusters * 2)
    , storedTile_(0)
    , binIndexMap_(matchDistribution, outputBinSize, skipEmptyBins)
    , flushThreads_(maxSavers)
    , binPathList_(buildBinPathList(binIndexMap_, matchDistribution, binDirectory,
                                    barcodeMetadataList, maxTileReads_, totalTiles, preSortBins, skipEmptyBins))
    , fragmentCollector_(binIndexMap_, maxTileClusters, flowcellLayoutList)
    , flushBuffer_(maxTileClusters, flowcellLayoutList)
    , threadDataFileBufCaches_(flushThreads_.size(),
                           FileBufCache(1, std::ios_base::out | std::ios_base::app | std::ios_base::binary))

{
    ISAAC_THREAD_CERR << "Resetting output files for " << binPathList_.size() << " bins" << std::endl;

    // assuming the last entry in the list contains the longest paths
    std::for_each(threadDataFileBufCaches_.begin(), threadDataFileBufCaches_.end(),
                  boost::bind(&FileBufCache::reservePathBuffers, _1,
                              binPathList_.back().getPathString().size()));

//    threadFragmentCollectors_.reserve(reserveClusters);

    ISAAC_THREAD_CERR << "Resetting output files done for " << binPathList_.size() << " bins" << std::endl;
}

void BufferingFragmentStorage::flushBin(
    const unsigned threadNumber,
    const unsigned binNumber, BinMetadata &binMetadata)
{
//    ISAAC_THREAD_CERR_DEV_TRACE((boost::format("BufferingFragmentStorage::flushBin: %d %s") % binNumber % binMetadata).str());

    ISAAC_ASSERT_MSG(binNumber == binMetadata.getIndex(), "Bin index mismatch");

//    ISAAC_THREAD_CERR_DEV_TRACE((boost::format("flushBuffer_.indexEnd() - currentBinIterator: %d") %
//        (flushBuffer_.indexEnd() - currentBinIterator)).str());

    const FragmentBuffer::IndexIterator binBegin = flushBuffer_.binBegin(binNumber, binIndexMap_);
    const FragmentBuffer::IndexIterator binEnd = flushBuffer_.binBegin(binNumber + 1, binIndexMap_);

    if (binBegin != binEnd && binBegin->initialized())
    {
        if (binMetadata.isEmpty())
        {
            // make sure file is empty first time we decide to put data in it.
            // boost::filesystem::remove for some stupid reason needs to allocate strings for this...
            unlink(binMetadata.getPath().c_str());
        }
        std::ostream osData(threadDataFileBufCaches_[threadNumber].get(binMetadata.getPath()));

        // store data sequentially in the bin file
        for(FragmentBuffer::IndexConstIterator currentBinIterator = binBegin;
            binEnd != currentBinIterator && currentBinIterator->initialized();
            ++currentBinIterator)
        {
            const FragmentBuffer::IndexRecord &recordStart = *currentBinIterator;

    //        ISAAC_THREAD_CERR << recordStart << std::endl;
    //        ISAAC_THREAD_CERR_DEV_TRACE((boost::format("%s") % recordStart).str());

            const io::FragmentHeader& header = recordStart.fragmentHeader();

            binMetadata.incrementDataSize(header.fStrandPosition_, header.getTotalLength());
    //        ISAAC_ASSERT_MSG(io::FragmentHeader::magicValue_ == header.magic_, "corrupt binary data in memory");
    //        ISAAC_ASSERT_MSG(header.getTotalLength() == header.totalLength_, "corrupt binary data in memory. fragment total length is bad");

            if (!osData.write(reinterpret_cast<const char *>(&header), header.getTotalLength())) {
                BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write into " + binMetadata.getPathString()));
            }

            if (!header.flags_.paired_)
            {
                binMetadata.incrementSeIdxElements(header.fStrandPosition_, 1, header.barcode_);
            }
            else if (header.flags_.reverse_ || header.flags_.unmapped_)
            {
                binMetadata.incrementRIdxElements(header.fStrandPosition_, 1, header.barcode_);
            }
            else
            {
                binMetadata.incrementFIdxElements(header.fStrandPosition_, 1, header.barcode_);
            }
            binMetadata.incrementGapCount(header.fStrandPosition_, header.gapCount_, header.barcode_);
            binMetadata.incrementCigarLength(header.fStrandPosition_, header.cigarLength_, header.barcode_);
        }
        // it is very important to flush these files after we're done. Although the FileBufWithReopen will do pubsync
        // itself, most likely the same file will be reopen on a different thread. this means, another FileBufWithReopen
        // will be already writing to it.
        osData.flush();
    }

//    ISAAC_THREAD_CERR << "Flushing bin done: " << binMetadata << std::endl;
}


void BufferingFragmentStorage::flushUnmappedBin(
    const unsigned threadNumber,
    const unsigned binNumber, BinMetadata &binMetadata)
{
//    ISAAC_THREAD_CERR_DEV_TRACE((boost::format("BufferingFragmentStorage::flushUnmappedBin: %d %s") % binNumber % binMetadata).str());

    ISAAC_ASSERT_MSG(binNumber == binMetadata.getIndex(), "Bin index mismatch");
    ISAAC_ASSERT_MSG(0 == binNumber, "Unaligned in index must be 0");

    const FragmentBuffer::IndexIterator binBegin = flushBuffer_.binBegin(binNumber, binIndexMap_);
    FragmentBuffer::IndexIterator binEnd = flushBuffer_.binBegin(binNumber + 1, binIndexMap_);
    if (binBegin != binEnd && binBegin->initialized())
    {

        if (binMetadata.isEmpty())
        {
            // make sure file is empty first time we decide to put data in it.
            // boost::filesystem::remove for some stupid reason needs to allocate strings for this...
            unlink(binMetadata.getPath().c_str());
        }
        std::ostream osData(threadDataFileBufCaches_[threadNumber].get(binMetadata.getPath()));

        FragmentBuffer::IndexIterator currentBinIterator = binBegin;
        unsigned long storedTileRead=0;
        for(;
            binEnd != currentBinIterator && currentBinIterator->initialized();
            ++currentBinIterator)
        {
            FragmentBuffer::IndexRecord &recordStart = *currentBinIterator;
            const io::FragmentAccessor& fragment = recordStart.fragment();
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
        // it is very important to flush these files after we're done. Although the FileBufWithReopen will do pubsync
        // itself, most likely the same file will be reopen on a different thread. this means, another FileBufWithReopen
        // will be already writing to it.
        osData.flush();
    }
//    ISAAC_THREAD_CERR << "Flushing bin done: " << binMetadata << std::endl;
}


void BufferingFragmentStorage::threadFlushBins(const unsigned threadNumber, unsigned &nextUnflushedBin)
{
    boost::lock_guard<boost::mutex> lock(binFlushMutex_);
    while(binPathList_.size() > nextUnflushedBin)
    {
        const unsigned ourBin = nextUnflushedBin++;

        {
            common::unlock_guard<boost::mutex > unlock(binFlushMutex_);

            BinMetadata &binMetadata = binPathList_.at(ourBin);
            if (0 == ourBin)
            {
                if (keepUnaligned_)
                {
                    flushUnmappedBin(threadNumber, ourBin, binMetadata);
                }
            }
            else
            {
                flushBin(threadNumber, ourBin, binMetadata);
            }
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

    ISAAC_ASSERT_MSG(binPathList_.size() == nextUnflushedBin, "Discrepancy between total bins and flushed bins count");

    flushBuffer_.clear();
    ++storedTile_;

    ISAAC_THREAD_CERR << "Flushing buffer done for " << nextUnflushedBin << " bins" << std::endl;
}

alignment::BinMetadataList BufferingFragmentStorage::buildBinPathList(
    const BinIndexMap &binIndexMap,
    const MatchDistribution &matchDistribution,
    const bfs::path &binDirectory,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const unsigned long maxTileReads,
    const unsigned long totalTiles,
    const bool preSortBins,
    const bool skipEmptyBins)
{
    ISAAC_THREAD_CERR << "maxTileClusters " << maxTileReads << "totalTiles " << totalTiles << std::endl;
    ISAAC_TRACE_STAT("before BufferingFragmentStorage::buildBinPathList");
    alignment::BinMetadataList binPathList;
    ISAAC_ASSERT_MSG(!binIndexMap.empty(), "Empty binIndexMap is illegal");
    ISAAC_ASSERT_MSG(!binIndexMap.back().empty(), "Empty binIndexMap entry is illegal" << binIndexMap);
    binPathList.reserve(1 + binIndexMap.back().back());
    size_t contigIndex = 0;
    BOOST_FOREACH(const std::vector<unsigned> &contigBins, binIndexMap)
    {
        ISAAC_ASSERT_MSG(!contigBins.empty(), "Unexpected empty contigBins");
        // matchDistribution contig 0 is the first contig
        // binIndexMap contig 0  is unaligned bin
        if (!skipEmptyBins || !contigIndex || !matchDistribution.isEmptyContig(contigIndex - 1))
        {
            for (unsigned i = contigBins.front(); contigBins.back() >= i; ++i)
            {
                ISAAC_ASSERT_MSG(binPathList.size() == i, "Basic sanity checking for bin numbering failed");
                using boost::format;
                const reference::ReferencePosition binStartPos = binIndexMap.getBinFirstPos(i);
                ISAAC_ASSERT_MSG(!i || binIndexMap.getBinIndex(binStartPos) == i, "BinIndexMap is broken");
                binPathList.push_back(
                    alignment::BinMetadata(
                        barcodeMetadataList.size(),
                        binPathList.size(),
                        binStartPos,
                        // bin zero has length of -1U as it contains unaligned records which are chunked by 32 bases of their sequence
                        i ? binIndexMap.getBinFirstInvalidPos(i) - binStartPos : maxTileReads * totalTiles,
                        // Pad file names well, so that we don't have to worry about them becoming of different length.
                        // This is important for memory reservation to be stable
                        binDirectory / (format("bin-%08d-%08d.dat") % contigIndex % i).str(),
                        /// Normally, aim to have 1024 or less chunks.
                        /// This will require about 4096*1024 (4 megabytes) of cache when pre-sorting bin during the loading in bam generator
                        preSortBins ? 1024 : 0));
            }
        }
        ++contigIndex;
    }
    return binPathList;
}

} //namespace matchSelector
} // namespace alignment
} // namespace isaac
