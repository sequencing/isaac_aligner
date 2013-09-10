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
 ** \file BinSorter.cpp
 **
 ** Reorders aligned data and stores results in bam file.
 ** 
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "bam/Bam.hh"

#include "build/BinSorter.hh"
#include "common/Memory.hh"

namespace isaac
{
namespace build
{

unsigned long BinSorter::serialize(
    boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams,
    boost::ptr_vector<bam::BamIndexPart> &bamIndexParts)
{
    ISAAC_THREAD_CERR << "Serializing records: " << getUniqueRecordsCount() <<  " of them for bin " << bin_ << std::endl;

    common::TimeSpec serTimeStart;
    ISAAC_ASSERT_MSG(-1 != clock_gettime(CLOCK_REALTIME, &serTimeStart), "clock_gettime failed, errno: " << errno << strerror(errno));

    if (isUnalignedBin())
    {
        unsigned long offset = 0;
        while(data_.size() != offset)
        {
            const io::FragmentAccessor &fragment = data_.getFragment(offset);
            bamSerializer_(fragment, bgzfStreams, bamIndexParts);
            offset += fragment.getTotalLength();
        }
    }
    else
    {
        const BaseType &v(*this);
        BOOST_FOREACH(const PackedFragmentBuffer::Index& idx, v)
        {
            bamSerializer_(idx, bgzfStreams, bamIndexParts, data_);
        }
    }

    BOOST_FOREACH(boost::iostreams::filtering_ostream &bgzfStream, bgzfStreams)
    {
        ISAAC_ASSERT_MSG(bgzfStream.strict_sync(), "Expecting the compressor to flush all the data");
    }

    common::TimeSpec serTimeEnd;
    ISAAC_ASSERT_MSG(-1 != clock_gettime(CLOCK_REALTIME, &serTimeEnd), "clock_gettime failed, errno: " << errno << strerror(errno));

    ISAAC_THREAD_CERR << "Serializing records done: " << getUniqueRecordsCount() <<  " of them for bin " << bin_ << " in " << common::tsdiff(serTimeStart, serTimeEnd) << "seconds." << std::endl;
    return size();
}

void verifyFragmentIntegrity(const io::FragmentAccessor &fragment)
{
/*
    if (fragment.getTotalLength() != fragment.totalLength_)
    {
        ISAAC_THREAD_CERR << pFragment << std::endl;
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Corrupt fragment (fragment total length is broken at %ld) read from %s") %
                (reinterpret_cast<const char *>(pFragment) - &data_.front()) %
                bin.getPathString()).str()));
    }
    if (io::FragmentAccessor::magicValue_ != fragment.magic_)
    {
        ISAAC_THREAD_CERR << pFragment << std::endl;
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Corrupt fragment (magic is broken) read from %s") % bin.getPathString()).str()));
    }
*/
}

void BinSorter::loadData()
{
    if(isUnalignedBin())
    {
        loadUnalignedData();
    }
    else
    {
        loadAlignedData();
    }
}

void BinSorter::loadUnalignedData()
{
    if(bin_.getDataSize())
    {
        ISAAC_THREAD_CERR << "Reading unaligned records from " << bin_ << std::endl;
        std::istream isData(fileBuf_.get(bin_.getPath()));
        if (!isData) {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open " + bin_.getPathString()));
        }

        if (!isData.seekg(bin_.getDataOffset(), std::ios_base::beg))
        {
            BOOST_THROW_EXCEPTION(common::IoException(
                errno, (boost::format("Failed to seek to position %d in %s") % bin_.getDataOffset() % bin_.getPathString()).str()));
        }
        if (!isData.read(&data_.front(), bin_.getDataSize())) {
            BOOST_THROW_EXCEPTION(common::IoException(
                errno, (boost::format("Failed to read %d bytes from %s") % bin_.getDataSize() % bin_.getPathString()).str()));
        }

/*
        unsigned count = 0;
        unsigned long offset = 0;
        while(data_.size() != offset)
        {
            const io::FragmentAccessor &fragment = data_.getFragment(offset);
            offset += fragment.getTotalLength();
            ++count;
        }
*/
        ISAAC_THREAD_CERR << "Reading unaligned records done from " << bin_ << std::endl;
    }
    else
    {
//        ISAAC_THREAD_CERR << "No unaligned records to read done from " << bin_ << std::endl;
    }
}


const io::FragmentAccessor &BinSorter::loadFragment(std::istream &isData, unsigned long &offset)
{
    io::FragmentHeader header;
    if (!isData.read(reinterpret_cast<char*>(&header), sizeof(header))) {
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Failed to read FragmentHeader bytes from %s") % bin_.getPathString()).str()));
    }

    const unsigned fragmentLength = header.getTotalLength();
    offset = dataDistribution_.addBytes(
        header.fStrandPosition_ - bin_.getBinStart(), fragmentLength);
//            ISAAC_THREAD_CERR << "offset:" << offset << " fragment: " << header << std::endl;
    io::FragmentAccessor &fragment = data_.getFragment(offset);
    io::FragmentHeader &headerRef = fragment;
    headerRef = header;
    if (!isData.read(reinterpret_cast<char*>(&fragment) + sizeof(header), fragmentLength - sizeof(header))) {
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Failed to read %d bytes from %s") % fragmentLength % bin_.getPathString()).str()));
    }

//    ISAAC_THREAD_CERR << "LOADED: " << fragment << std::endl;

    return fragment;
}

void BinSorter::loadAlignedData()
{
    if(bin_.getDataSize())
    {
        ISAAC_THREAD_CERR << "Reading alignment records from " << bin_ << std::endl;
        unsigned long dataSize = 0;
        // summarize chunk sizes to get offsets
        dataDistribution_.tallyOffsets();
        std::istream isData(fileBuf_.get(bin_.getPath()));
        if (!isData) {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open " + bin_.getPathString()));
        }
        if (!isData.seekg(bin_.getDataOffset()))
        {
            BOOST_THROW_EXCEPTION(common::IoException(
                errno, (boost::format("Failed to seek to position %d in %s") % bin_.getDataOffset() % bin_.getPathString()).str()));
        }

        rIdxFileContent_.clear();
        fIdxFileContent_.clear();
        seIdxFileContent_.clear();
        unsigned long binOffset = 0;

        while(isData && dataSize != bin_.getDataSize())
        {
            unsigned long offset = 0;
            const io::FragmentAccessor &fragment = loadFragment(isData, offset);
            const unsigned fragmentLength = fragment.getTotalLength();
            dataSize += fragmentLength;

            verifyFragmentIntegrity(fragment);

            if (!fragment.flags_.paired_)
            {
                SeFragmentIndex seIdx(fragment.fStrandPosition_);
                seIdx.dataOffset_ = offset;
                seIdxFileContent_.push_back(seIdx);
                binOffset += fragmentLength;
            }
            else
            {
                binOffset += fragmentLength;
                unsigned long mateOffset = offset;
                if (bin_.coversPosition(fragment.mateFStrandPosition_))
                {
                    const io::FragmentAccessor &mateFragment = loadFragment(isData, mateOffset);
                    ISAAC_ASSERT_MSG(mateFragment.clusterId_ == fragment.clusterId_, "mateFragment.clusterId_ != fragment.clusterId_");
                    ISAAC_ASSERT_MSG(mateFragment.flags_.unmapped_ == fragment.flags_.mateUnmapped_, "mateFragment.flags_.unmapped_ != fragment.flags_.mateUnmapped_");
                    ISAAC_ASSERT_MSG(mateFragment.flags_.reverse_ == fragment.flags_.mateReverse_,
                                     "mateFragment.flags_.reverse_ != fragment.flags_.mateReverse_" << fragment << " " << mateFragment);

                    const unsigned mateFragmentLength = mateFragment.getTotalLength();
                    dataSize += mateFragmentLength;

                    verifyFragmentIntegrity(mateFragment);

                    if (mateFragment.flags_.reverse_ || mateFragment.flags_.unmapped_)
                    {
                        RStrandOrShadowFragmentIndex rsMateIdx(
                            mateFragment.fStrandPosition_, // shadows are stored at the position of their singletons,
                            io::FragmentIndexAnchor(mateFragment),
                            FragmentIndexMate(
                                mateFragment.flags_.mateUnmapped_, mateFragment.flags_.mateReverse_, mateFragment.mateStorageBin_,
                                mateFragment.mateAnchor_),
                                mateFragment.duplicateClusterRank_);

                        rsMateIdx.dataOffset_ = mateOffset;
                        rsMateIdx.mateDataOffset_ = offset;
                        rIdxFileContent_.push_back(rsMateIdx);
                    }
                    else
                    {
                        FStrandFragmentIndex fMateIdx(
                            mateFragment.fStrandPosition_,
                            FragmentIndexMate(
                                mateFragment.flags_.mateUnmapped_, mateFragment.flags_.mateReverse_, mateFragment.mateStorageBin_,
                                mateFragment.mateAnchor_),
                                mateFragment.duplicateClusterRank_);

                        fMateIdx.dataOffset_ = mateOffset;
                        fMateIdx.mateDataOffset_ = offset;
                        fIdxFileContent_.push_back(fMateIdx);
                    }

                    binOffset += mateFragmentLength;
                }

                if (fragment.flags_.reverse_ || fragment.flags_.unmapped_)
                {
                    RStrandOrShadowFragmentIndex rsIdx(
                        fragment.fStrandPosition_, // shadows are stored at the position of their singletons,
                        io::FragmentIndexAnchor(fragment),
                        FragmentIndexMate(
                            fragment.flags_.mateUnmapped_, fragment.flags_.mateReverse_, fragment.mateStorageBin_,
                            fragment.mateAnchor_),
                        fragment.duplicateClusterRank_);

                    rsIdx.dataOffset_ = offset;
                    rsIdx.mateDataOffset_ = mateOffset;
                    rIdxFileContent_.push_back(rsIdx);
                }
                else
                {
                    FStrandFragmentIndex fIdx(
                        fragment.fStrandPosition_,
                        FragmentIndexMate(
                            fragment.flags_.mateUnmapped_, fragment.flags_.mateReverse_, fragment.mateStorageBin_,
                            fragment.mateAnchor_),
                        fragment.duplicateClusterRank_);

                    fIdx.dataOffset_ = offset;
                    fIdx.mateDataOffset_ = mateOffset;
                    fIdxFileContent_.push_back(fIdx);
                }
            }
        }
        ISAAC_THREAD_CERR << "Reading alignment records done from " << bin_ << std::endl;
    }
}


void BinSorter::resolveDuplicates(
    BuildStats &buildStats)
{
    NotAFilter().filterInput(data_, seIdxFileContent_.begin(), seIdxFileContent_.end(), buildStats, binStatsIndex_, std::back_inserter<BaseType>(*this));
    if (keepDuplicates_ && !markDuplicates_)
    {
        NotAFilter().filterInput(data_, rIdxFileContent_.begin(), rIdxFileContent_.end(), buildStats, binStatsIndex_, std::back_inserter<BaseType>(*this));
        NotAFilter().filterInput(data_, fIdxFileContent_.begin(), fIdxFileContent_.end(), buildStats, binStatsIndex_, std::back_inserter<BaseType>(*this));
    }
    else
    {
        if (singleLibrarySamples_)
        {
            DuplicatePairEndFilter(keepDuplicates_).filterInput(
                RSDuplicateFilter<true>(barcodeBamMapping_.getSampleIndexMap()),
                data_, rIdxFileContent_.begin(), rIdxFileContent_.end(),
                buildStats, binStatsIndex_, std::back_inserter<BaseType>(*this));
            DuplicatePairEndFilter(keepDuplicates_).filterInput(
                FDuplicateFilter<true>(barcodeBamMapping_.getSampleIndexMap()),
                data_, fIdxFileContent_.begin(), fIdxFileContent_.end(),
                buildStats, binStatsIndex_, std::back_inserter<BaseType>(*this));
        }
        else
        {
            DuplicatePairEndFilter(keepDuplicates_).filterInput(
                RSDuplicateFilter<false>(barcodeBamMapping_.getSampleIndexMap()),
                data_, rIdxFileContent_.begin(), rIdxFileContent_.end(),
                buildStats, binStatsIndex_, std::back_inserter<BaseType>(*this));
            DuplicatePairEndFilter(keepDuplicates_).filterInput(
                FDuplicateFilter<false>(barcodeBamMapping_.getSampleIndexMap()),
                data_, fIdxFileContent_.begin(), fIdxFileContent_.end(),
                buildStats, binStatsIndex_, std::back_inserter<BaseType>(*this));
        }
    }
}

inline size_t getTotalGapsCount(const std::vector<RealignerGaps> &allGaps)
{
    return std::accumulate(allGaps.begin(), allGaps.end(), 0,
                           boost::bind(std::plus<size_t>(), _1,
                                       boost::bind(&RealignerGaps::getGapsCount, _2)));
}

unsigned BinSorter::getGapGroupsCount() const
{
    switch (realignGaps_)
    {
    case REALIGN_SAMPLE:
        return barcodeBamMapping_.getTotalSamples();
    case REALIGN_PROJECT:
        return barcodeBamMapping_.getMaxProjectIndex() + 1;
    case REALIGN_ALL:
        return 1;
    case REALIGN_NONE:
        return 0;
    default:
        ISAAC_ASSERT_MSG(false, "Unknown gap realignment mode: " << realignGaps_);
        break;
    }
    return -1;
}

inline unsigned BinSorter::getGapGroupIndex(const unsigned barcode) const
{
    switch (realignGaps_)
    {
    case REALIGN_SAMPLE:
        return barcodeBamMapping_.getSampleIndex(barcode);
    case REALIGN_PROJECT:
        return barcodeBamMapping_.getProjectIndex(barcode);
    case REALIGN_ALL:
        return 0;
    default:
        ISAAC_ASSERT_MSG(false, "Unknown gap realignment mode: " << realignGaps_);
        break;
    }
    return -1;
}

void BinSorter::reserveGaps(
    const alignment::BinMetadata& bin,
    const flowcell::BarcodeMetadataList &barcodeMetadataList)
{
    std::vector<size_t> gapsByGroup(getGapGroupsCount(), 0);
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        gapsByGroup.at(getGapGroupIndex(barcode.getIndex())) +=
            bin.getBarcodeGapCount(barcode.getIndex());
    }
    unsigned gapGroupId = 0;
    BOOST_FOREACH(const size_t gaps, gapsByGroup)
    {
        realignerGaps_.at(gapGroupId++).reserve(gaps);
    }
}

void BinSorter::collectGaps()
{
    for(std::vector<char>::const_iterator p = data_.begin(); p != data_.end();)
    {
        const io::FragmentAccessor &fragment = *reinterpret_cast<const io::FragmentAccessor *>(&*p);

        if (fragment.gapCount_)
        {
            const unsigned gapGroupIndex = getGapGroupIndex(fragment.barcode_);
            realignerGaps_.at(gapGroupIndex).addGaps(fragment.fStrandPosition_, fragment.cigarBegin(), fragment.cigarEnd());
        }
        p += fragment.getTotalLength();
    }
    ISAAC_THREAD_CERR << "Finalizing " << getTotalGapsCount(realignerGaps_) << " gaps."  << std::endl;

    std::for_each(realignerGaps_.begin(), realignerGaps_.end(), boost::bind(&RealignerGaps::finalizeGaps, _1));
}

void BinSorter::realignGaps()
{
    ISAAC_THREAD_CERR << "Realigning against " << getTotalGapsCount(realignerGaps_) << " unique gaps. " << bin_ << std::endl;
    BOOST_FOREACH(PackedFragmentBuffer::Index &index, std::make_pair(indexBegin(), indexEnd()))
    {
        io::FragmentAccessor &fragment = data_.getFragment(index);
        const unsigned gapGroupIndex = getGapGroupIndex(fragment.barcode_);

        gapRealigner_.realign(realignerGaps_.at(gapGroupIndex), bin_.getBinStart(), bin_.getBinEnd(), index, fragment, data_);
    }
    ISAAC_THREAD_CERR << "Realigning gaps done" << std::endl;
}

} // namespace build
} // namespace isaac
