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
 ** \file BgzfReader.cpp
 **
 ** Component to read bgzf blocks.
 **
 ** \author Roman Petrovski
 **/

#include <boost/format.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>


#include "bgzf/BgzfReader.hh"

namespace isaac
{
namespace bgzf
{

inline void validateHeader(const bgzf::Header &header)
{
    ISAAC_ASSERT_MSG(header.ID1 == 31U, " got " << unsigned(header.ID1));
    ISAAC_ASSERT_MSG(header.ID2 == 139U, " got " << unsigned(header.ID2));
    ISAAC_ASSERT_MSG(header.CM == 8U, " got " << header.CM);
    ISAAC_ASSERT_MSG(header.xfield.XLEN[0] == 6U, " got " << unsigned(header.xfield.XLEN[0]));
    ISAAC_ASSERT_MSG(header.xfield.XLEN[1] == 0U, " got " << unsigned(header.xfield.XLEN[1]));
    ISAAC_ASSERT_MSG(header.xfield.SI1 == 66U, " got " << unsigned(header.xfield.SI1));
    ISAAC_ASSERT_MSG(header.xfield.SI2 == 67U, " got " << unsigned(header.xfield.SI2));
}

unsigned BgzfReader::readNextBlock(std::istream &is)
{
    compressedBlockBuffer_.clear();

    unsigned ret = 0;
    for (unsigned i = 0; i < BLOCKS_AT_ONCE; ++i)
    {
        compressedBlockBuffer_.resize(compressedBlockBuffer_.size() + sizeof(bgzf::Header));
        is.read(&compressedBlockBuffer_.front() + compressedBlockBuffer_.size() - sizeof(bgzf::Header), sizeof(bgzf::Header));
        if (is.eof())
        {
            break;
        }
        if (!is)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read bam file header: %s") %
                strerror(errno)).str()));
        }
        const bgzf::Header &header = *reinterpret_cast<bgzf::Header *>(
            &compressedBlockBuffer_.front() + compressedBlockBuffer_.size() - sizeof(bgzf::Header));
        validateHeader(header);

        compressedBlockBuffer_.resize(compressedBlockBuffer_.size() + header.getCDATASize() + sizeof(bgzf::Footer));
        if (!is.read(&compressedBlockBuffer_.front() + compressedBlockBuffer_.size() - header.getCDATASize() - sizeof(bgzf::Footer), header.getCDATASize() + sizeof(bgzf::Footer)))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read %d bytes of bam CDATA : %s") %
                header.getCDATASize() % strerror(errno)).str()));
        }
    //    ISAAC_THREAD_CERR << "header.getCDATASize() " << header.getCDATASize() << " bytes" << std::endl;

        if (header.getCDATASize())
        {
            const bgzf::Footer &footer = *reinterpret_cast<bgzf::Footer*>(
                &compressedBlockBuffer_.front() + compressedBlockBuffer_.size() - sizeof(bgzf::Footer));
    //        ISAAC_ASSERT_MSG(footer.getISIZE(), "Unexpected ISIZE 0 for CDATA size " << header.getCDATASize());
            ret += footer.getISIZE();
        }
    }

    return ret;
}

void BgzfReader::uncompressCurrentBlock(char* p, const std::size_t size)
{
    reset();
    strm_.next_out = reinterpret_cast<Bytef *>(p);
    strm_.avail_out = size;

    strm_.next_in = reinterpret_cast<Bytef *>(&compressedBlockBuffer_.front());
    strm_.avail_in = compressedBlockBuffer_.size();
    int err = inflate(&strm_, Z_SYNC_FLUSH);
    if (Z_OK != err && Z_STREAM_END != err)
    {
        BOOST_THROW_EXCEPTION(BgzfInflateException(err, strm_));
    }
    const std::size_t decompressedBytes = size - strm_.avail_out;

    boost::iostreams::basic_array_source<char> src(&compressedBlockBuffer_.front(), compressedBlockBuffer_.size());

    if (decompressedBytes != size)
    {
        BOOST_THROW_EXCEPTION(common::IoException(EINVAL, (boost::format("Unexpected number of BGZF bytes uncompressed. "
            "Expected %d, uncompressed %d") % size % decompressedBytes).str()));
    }
}

void ParallelBgzfReader::open(const boost::filesystem::path &bamPath)
{
    fileBuffer_.reopen(bamPath.c_str(), io::FileBufWithReopen::SequentialOnce);
    if (!fileBuffer_.is_open())
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to open bam file: %s") % bamPath).str()));
    }
    pendingBlockSize_ = 0;
    is_.rdbuf(&fileBuffer_);
    ISAAC_THREAD_CERR << "Opened bam stream on " << bamPath << std::endl;
}

void ParallelBgzfReader::waitForLoadSlot(boost::unique_lock<boost::mutex> &lock)
{
    while (!loadSlotAvailable_)
    {
        stateChangedCondition_.wait(lock);
    }
    loadSlotAvailable_ = false;
}

void ParallelBgzfReader::releaseLoadSlot()
{
    ISAAC_ASSERT_MSG(!loadSlotAvailable_, "Invalid attempt to release a load slot that is not locked");
    loadSlotAvailable_ = true;
    stateChangedCondition_.notify_all();
}

void ParallelBgzfReader::waitForComputeSlot(boost::unique_lock<boost::mutex> &lock)
{
    while (!computeSlotsAvailable_)
    {
        stateChangedCondition_.wait(lock);
    }
    --computeSlotsAvailable_;
}

void ParallelBgzfReader::releaseComputeSlot()
{
    ++computeSlotsAvailable_;
    stateChangedCondition_.notify_all();
}

void ParallelBgzfReader::readMoreDataParallel(const unsigned threadNumber, std::vector<char> &buffer)
{
    boost::unique_lock<boost::mutex> lock(stateMutex_);
    while (true)
    {
        unsigned long &ourThreadOffset = threadOffsets_.at(threadNumber);
        unsigned long &ourThreadBlockSize = threadBlockSizes_.at(threadNumber);

        // first get rid of any available data that is pending delivery
        if (-1UL != ourThreadOffset)
        {
            ISAAC_ASSERT_MSG(buffer.size() >= ourThreadOffset + ourThreadBlockSize,
                             "Result buffer is unexpectedly insufficient to fit the uncompressed block."
                             " ourThreadOffset: " << ourThreadOffset <<
                             " ourThreadBlockSize" << ourThreadBlockSize <<
                             " buffer.size()" << buffer.size());
            waitForComputeSlot(lock);
            ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&ParallelBgzfReader::releaseComputeSlot, this))
            {
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                readers_.at(threadNumber).uncompressCurrentBlock(&buffer.front() + ourThreadOffset, ourThreadBlockSize);
            }
            ourThreadOffset = -1UL;
        }

        if (isEof())
        {
            ISAAC_THREAD_CERR << "Thread " << threadNumber << " terminating due to eof" << " buffer.size()" << buffer.size() << std::endl;
            break;
        }

        // Load next bgzf block
        BgzfReader &ourThreadReader = readers_.at(threadNumber);
        waitForLoadSlot(lock);
        ourThreadBlockSize = 0;
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&ParallelBgzfReader::releaseLoadSlot, this))
        {
            if (nextUncompressedOffset_ < buffer.capacity())
            {
                {
                    common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                    while (!isEof() && !ourThreadBlockSize)
                    {
                        ourThreadBlockSize = ourThreadReader.readNextBlock(is_);
                    }
                }
                if (!ourThreadBlockSize)
                {
                    ISAAC_ASSERT_MSG(isEof(), "Unexpectedly stopped reading bam before the end of the input stream");
                    ISAAC_THREAD_CERR << "Thread " << threadNumber << " reached eof while reading empty blocks" << std::endl;
                    break;
                }
                else if (isEof())
                {
                    ISAAC_THREAD_CERR << "Thread " << threadNumber << " reached eof" << std::endl;
                }
//                ISAAC_THREAD_CERR << "Thread " << threadNumber << " read " << ourThreadBlockSize << " bytes" << std::endl;
                ourThreadOffset = nextUncompressedOffset_;
                nextUncompressedOffset_ += ourThreadBlockSize;

                if (nextUncompressedOffset_ <= buffer.capacity())
                {
                    buffer.resize(nextUncompressedOffset_);
                    pendingBlockSize_ = 0;
                }
                else
                {
                    ourThreadOffset = 0;
                    pendingBlockSize_ = ourThreadBlockSize;
//                    ISAAC_THREAD_CERR << "Thread " << threadNumber << " has no room to place " << ourThreadBlockSize << " bytes" << std::endl;
                    break;
                }
            }
            else
            {
//                ISAAC_THREAD_CERR << "Thread " << threadNumber <<
//                    " nextUncompressedOffset_ " << nextUncompressedOffset_ <<
//                    " buffer.capacity()" << buffer.capacity() << std::endl;
                break;
            }
        }
    }

}

bool ParallelBgzfReader::readMoreData(std::vector<char> &buffer)
{
    std::vector<unsigned long>::iterator pendingBlockOffsetIt =
        std::find(threadOffsets_.begin(), threadOffsets_.end(), 0UL);
    if (threadOffsets_.end() != pendingBlockOffsetIt)
    {
        // preserve existing buffer.size() as there is some data the client wants to merge with what we're about
        // to decompress
        *pendingBlockOffsetIt = buffer.size();
    }
    nextUncompressedOffset_ = pendingBlockSize_ + buffer.size();
    pendingBlockSize_ = 0;
    buffer.resize(nextUncompressedOffset_);
//    ISAAC_THREAD_CERR << "ParallelBgzfReader::readMoreData "
//        " pendingBlockSize_" << pendingBlockSize_ <<
//        " buffer.size()" << buffer.size() << std::endl;
    const std::size_t oldSize = buffer.size();
    threads_.execute(boost::bind(&ParallelBgzfReader::readMoreDataParallel, this, _1, boost::ref(buffer)), coresMax_);
    return threadOffsets_.end() != pendingBlockOffsetIt || oldSize != buffer.size();
}


} // namespace bgzf
} // namespace isaac
