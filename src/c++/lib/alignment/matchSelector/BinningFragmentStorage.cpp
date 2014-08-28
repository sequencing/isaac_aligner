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
 ** \file BinningFragmentStorage.cpp
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
#include "alignment/matchSelector/BinningFragmentStorage.hh"
#include "io/Fragment.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

BinningFragmentStorage::BinningFragmentStorage(
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
    const unsigned long totalTiles)
    : keepUnaligned_(keepUnaligned)
    , maxTileReads_(maxTileClusters * READS_MAX)
    , binIndexMap_(matchDistribution, outputBinSize, false)
    , binPathList_(buildBinPathList(binIndexMap_, matchDistribution.getBinSize(), binDirectory,
                                    barcodeMetadataList, maxTileReads_, totalTiles, preSortBins))
    , binFiles_(binPathList_.size())

{
    ISAAC_THREAD_CERR << "Resetting output files for " << binPathList_.size() << " bins" << std::endl;

    BOOST_FOREACH(const BinMetadata &binMetadata, binPathList_)
    {
        binFiles_.push_back(new std::ofstream(binMetadata.getPath().c_str(), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary));

        if (!binFiles_.back())
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open bin file " + binMetadata.getPathString()));
        }
    }

    ISAAC_THREAD_CERR << "Resetting output files done for " << binPathList_.size() << " bins" << std::endl;
}

alignment::BinMetadataList BinningFragmentStorage::buildBinPathList(
    const BinIndexMap &binIndexMap,
    const unsigned long outputBinSize,
    const bfs::path &binDirectory,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const unsigned long maxTileReads,
    const unsigned long totalTiles,
    const bool preSortBins)
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
            // Pad file names well, so that we don't have to worry about them becoming of different length.
            // This is important for memory reservation to be stable
            const reference::ReferencePosition binStartPos = binIndexMap.getBinFirstPos(i);
            binPathList.push_back(
                alignment::BinMetadata(
                    barcodeMetadataList.size(),
                    binPathList.size(),
                    binStartPos,
                    // bin zero has length of -1U as it contains unaligned records which are chunked by 32 bases of their sequence
                    i ? binIndexMap.getBinFirstInvalidPos(i) - binStartPos : maxTileReads * totalTiles,
                    binDirectory / (format("bin-%04d-%04d.dat") % contigIndex % i).str(),
                    // Normally, aim to have 1024 or less chunks.
                    // This will require about 4096*1024 (4 megabytes) of cache when pre-sorting bin during the loading in bam generator
                    preSortBins ? 1024 : 0));
        }
        ++contigIndex;
    }
    return binPathList;
}

template <typename BufferT>
void storeBclAndCigar(
    const alignment::FragmentMetadata & fragment,
    BufferT &buffer)
{
    // copy the bcl data (reverse-complement the sequence if the fragment is reverse-aligned)
    std::vector<char>::const_iterator bclData = fragment.getBclData();
    if (fragment.isReverse())
    {
        std::transform(std::reverse_iterator<std::vector<char>::const_iterator>(bclData + fragment.getReadLength()),
                                      std::reverse_iterator<std::vector<char>::const_iterator>(bclData),
                                      std::back_inserter(buffer), oligo::getReverseBcl);
    }
    else
    {
        std::copy(bclData, bclData + fragment.getReadLength(), std::back_inserter(buffer));

    }

    if (fragment.isAligned())
    {
        const alignment::Cigar::const_iterator cigarBegin = fragment.cigarBuffer->begin() + fragment.cigarOffset;
        const alignment::Cigar::const_iterator cigarEnd = cigarBegin + fragment.cigarLength;
        BOOST_FOREACH(const unsigned cig, std::make_pair(cigarBegin, cigarEnd))
        {
            buffer.push_back(cig);
            buffer.push_back(cig >> 8);
            buffer.push_back(cig >> 16);
            buffer.push_back(cig >> 24);
        }
    }
}

static const unsigned FRAGMENT_BYTES_MAX = 10*1024;

template <typename BufferT>
unsigned BinningFragmentStorage::packFragment(
    const alignment::BamTemplate &bamTemplate,
    const unsigned fragmentIndex,
    const unsigned barcodeIdx,
    BufferT &buffer)
{
    ISAAC_ASSERT_MSG(READS_MAX >= bamTemplate.getFragmentCount(), "Expected paired or single-ended data");

    const alignment::FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(fragmentIndex);

    unsigned storageBin = -1;

    const alignment::FragmentMetadata &mate = bamTemplate.getMateFragmentMetadata(fragment);
    if (mate.isNoMatch() && fragment.isNoMatch())
    {
        storageBin = 0;
        const io::FragmentHeader header(bamTemplate, fragment, mate, barcodeIdx, 0);
        std::copy(header.bytesBegin(), header.bytesEnd(), std::back_inserter(buffer));
    }
    else
    {
        const unsigned mateStorageBin = binIndexMap_.getBinIndex(mate.getFStrandReferencePosition());
        const io::FragmentHeader header(bamTemplate, fragment, mate, barcodeIdx,
                                                          mateStorageBin);

        storageBin = binIndexMap_.getBinIndex(header.fStrandPosition_);

        std::copy(header.bytesBegin(), header.bytesEnd(), std::back_inserter(buffer));
    }


    storeBclAndCigar(fragment, buffer);
    ISAAC_ASSERT_MSG(buffer.size() == reinterpret_cast<io::FragmentHeader&>(buffer.front()).getTotalLength(),
                     "buffer.size()=" << buffer.size() << " " << reinterpret_cast<io::FragmentHeader&>(buffer.front()));
    return storageBin;
}

template <typename BufferT>
void BinningFragmentStorage::storeFragment(
    const BufferT &buffer,
    const unsigned storageBin)
{
    const io::FragmentHeader &header = reinterpret_cast<const io::FragmentHeader &>(buffer.front());
    ISAAC_ASSERT_MSG(buffer.size() == header.getTotalLength(), "buffer.size()=" << buffer.size() << " " << header);

    BinMetadata &binMetadata = binPathList_.at(storageBin);
    if (!header.isAligned() && !header.isMateAligned())
    {
        const unsigned long globalReadId = maxTileReads_ * header.tile_ +
            header.clusterId_ * 2 + header.flags_.secondRead_;
        binMetadata.incrementDataSize(globalReadId, header.getTotalLength());
        binMetadata.incrementNmElements(globalReadId, 1, header.barcode_);
    }
    else
    {
        binMetadata.incrementDataSize(header.fStrandPosition_, header.getTotalLength());
        if (header.flags_.reverse_ || header.flags_.unmapped_)
        {
            binMetadata.incrementRIdxElements(header.fStrandPosition_, 1, header.barcode_);
        }
        else
        {
            binMetadata.incrementFIdxElements(header.fStrandPosition_, 1, header.barcode_);
        }
        binMetadata.incrementGapCount(header.fStrandPosition_, header.gapCount_, header.barcode_);
        binMetadata.incrementCigarLength(header.fStrandPosition_, header.cigarLength_, header.barcode_);
        ISAAC_ASSERT_MSG(binMetadata.getBinStart().getContigId() == header.fStrandPosition_.getContigId(), "tada: " << binMetadata << header);
    }

    std::ostream &osData = binFiles_.at(storageBin);
    if (!osData.write(&buffer.front(), buffer.size())) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write into " + binMetadata.getPathString()));
    }
}

void BinningFragmentStorage::add(const BamTemplate &bamTemplate, const unsigned barcodeIdx)
{
    if (2 == bamTemplate.getFragmentCount())
    {
        storePaired(bamTemplate, barcodeIdx);
    }
    else
    {
        storeSingle(bamTemplate, barcodeIdx);
    }
}

void BinningFragmentStorage::storeSingle(const BamTemplate &bamTemplate, const unsigned barcodeIdx)
{
    common::FiniteCapacityVector<char, sizeof(io::FragmentHeader) + FRAGMENT_BYTES_MAX> buffer;
    const alignment::FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(0);
    const unsigned storageBin = fragment.isNoMatch() ? 0 : binIndexMap_.getBinIndex(fragment.getFStrandReferencePosition());

    const io::FragmentHeader header(bamTemplate, fragment, barcodeIdx);
    std::copy(header.bytesBegin(), header.bytesEnd(), std::back_inserter(buffer));

    boost::unique_lock<boost::mutex> lock(binMutex_[storageBin % binMutex_.size()]);

    if (fragment.isNoMatch())
    {
        BinMetadata &binMetadata = binPathList_.at(storageBin);
        const unsigned long globalReadId = maxTileReads_ * fragment.getCluster().getTile() +
            fragment.getCluster().getId();
        binMetadata.incrementDataSize(globalReadId, header.getTotalLength());
        binMetadata.incrementNmElements(globalReadId, 1, header.barcode_);
    }
    else
    {
        BinMetadata &binMetadata = binPathList_.at(storageBin);
        binMetadata.incrementDataSize(header.fStrandPosition_, header.getTotalLength());
        binMetadata.incrementSeIdxElements(header.fStrandPosition_, 1, header.barcode_);
        binMetadata.incrementGapCount(header.fStrandPosition_, header.gapCount_, header.barcode_);
        binMetadata.incrementCigarLength(header.fStrandPosition_, header.cigarLength_, header.barcode_);
    }

    storeBclAndCigar(fragment, buffer);
    ISAAC_ASSERT_MSG(buffer.size() == reinterpret_cast<io::FragmentHeader&>(buffer.front()).getTotalLength(),
                     "buffer.size()=" << buffer.size() << " " << reinterpret_cast<io::FragmentHeader&>(buffer.front()));

    std::ostream &osData = binFiles_.at(storageBin);
    if (!osData.write(&buffer.front(), buffer.size())) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write into " + binPathList_.at(storageBin).getPathString()));
    }
}

void BinningFragmentStorage::storePaired(const BamTemplate &bamTemplate, const unsigned barcodeIdx)
{
    common::FiniteCapacityVector<char, sizeof(io::FragmentHeader) + FRAGMENT_BYTES_MAX> bufferA;
    common::FiniteCapacityVector<char, sizeof(io::FragmentHeader) + FRAGMENT_BYTES_MAX> bufferB;

    const unsigned binA = packFragment(bamTemplate, 0, barcodeIdx, bufferA);
    const unsigned binB = packFragment(bamTemplate, 1, barcodeIdx, bufferB);

    {
        boost::unique_lock<boost::mutex> lock(binMutex_[binA % binMutex_.size()]);
        storeFragment(bufferA, binA);
        if (binB == binA)
        {
            // BinSorter requires that records follow each other if they are in the same bin.
            // ensure lock does not get released in between the two!
            storeFragment(bufferB, binB);
        }
    }

    if (binB != binA)
    {
        boost::unique_lock<boost::mutex> lock(binMutex_[binB % binMutex_.size()]);
        storeFragment(bufferB, binB);
    }
}

} //namespace matchSelector
} // namespace alignment
} // namespace isaac
