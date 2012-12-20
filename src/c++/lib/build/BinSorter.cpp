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
        ISAAC_THREAD_CERR << "No unaligned records to read done from " << bin_ << std::endl;
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

    return fragment;
}

void BinSorter::loadAlignedData()
{
    ISAAC_THREAD_CERR << "Reading alignment records from " << bin_ << std::endl;
    if(bin_.getDataSize())
    {
        unsigned long dataSize = 0;
        dataDistribution_ = bin_.getDataDistribution();
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

        std::vector<io::RStrandOrShadowFragmentIndex>::iterator rIdxIterator= rIdxFileContent_.begin();
        std::vector<io::FStrandFragmentIndex>::iterator fIdxIterator= fIdxFileContent_.begin();
        std::vector<io::SeFragmentIndex>::iterator seIdxIterator= seIdxFileContent_.begin();
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
                ISAAC_ASSERT_MSG(seIdxIterator->dataOffset_ == binOffset, "Unexpected offset in loaded se index" << fragment);
                seIdxIterator->dataOffset_ = offset;
                seIdxIterator++;
                binOffset += fragmentLength;
            }
            else
            {
                io::PairEndIndex *idxPointer = 0;
                if (fragment.flags_.reverse_ || fragment.flags_.unmapped_)
                {
                    ISAAC_ASSERT_MSG(rIdxIterator->dataOffset_ == binOffset, "Unexpected offset in loaded rs index" << fragment);
                    idxPointer = &*rIdxIterator++;
                }
                else
                {
                    ISAAC_ASSERT_MSG(fIdxIterator->dataOffset_ == binOffset, "Unexpected offset in loaded fw index" << fragment);
                    idxPointer = &*fIdxIterator++;
                }
                if (!idxPointer)
                {
                    ISAAC_ASSERT_MSG(false, "None of loaded index offsets matched the file offset for fragment: " << fragment << *fIdxIterator << *rIdxIterator << " binOffset:" << binOffset);
                }
                binOffset += fragmentLength;

                if (idxPointer->mateDataOffset_ != idxPointer->dataOffset_)
                //if (idxPointer->mate_.info_.fields_.storageBin_ == bin_.getIndex()) NOTE! the storageBin_ index is most likely not valid at this point
                    // due to the unaligned bin splitting
                {
                    ISAAC_ASSERT_MSG(binOffset == idxPointer->mateDataOffset_,
                                     "Mate fragments must be adjacent in the same bin. Expected " << binOffset << " actual:" << idxPointer->mateDataOffset_);

                    unsigned long mateOffset = 0;
                    const io::FragmentAccessor &mateFragment = loadFragment(isData, mateOffset);
                    const unsigned mateFragmentLength = mateFragment.getTotalLength();
                    dataSize += mateFragmentLength;

                    verifyFragmentIntegrity(mateFragment);

                    io::PairEndIndex *mateIdxPointer = 0;
                    if (mateFragment.flags_.reverse_ || mateFragment.flags_.unmapped_)
                    {
                        ISAAC_ASSERT_MSG(rIdxIterator->dataOffset_ == binOffset, "Unexpected offset in loaded rs mate index" << mateFragment);
                        mateIdxPointer = &*rIdxIterator++;
                    }
                    else
                    {
                        ISAAC_ASSERT_MSG(fIdxIterator->dataOffset_ == binOffset, "Unexpected offset in loaded fw mate index" << mateFragment);
                        mateIdxPointer = &*fIdxIterator++;
                    }

                    if (!mateIdxPointer)
                    {
                        ISAAC_ASSERT_MSG(false, "None of loaded index offsets matched the mate file offset for fragment: " << fragment << "-" << mateFragment);
                    }

                    mateIdxPointer->dataOffset_ = mateOffset;
                    mateIdxPointer->mateDataOffset_ = offset;
                    idxPointer->mateDataOffset_ = mateOffset;
                    binOffset += mateFragmentLength;
                }
                else
                {
                    idxPointer->mateDataOffset_ = offset;
                }

                idxPointer->dataOffset_ = offset;
            }
        }
    }
    ISAAC_THREAD_CERR << "Reading alignment records done from " << bin_ << std::endl;
}

} // namespace build
} // namespace isaac
