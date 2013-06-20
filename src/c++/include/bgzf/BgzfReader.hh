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
 ** \file BgzfReader.hh
 **
 ** Component to read bgzf blocks.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BGZF_BGZF_READER_HH
#define iSAAC_BGZF_BGZF_READER_HH

#include <zlib.h>

#include <boost/filesystem.hpp>

#include "bgzf/Bgzf.hh"
#include "common/FiniteCapacityVector.hh"
#include "common/Threads.hpp"
#include "io/FileBufWithReopen.hh"

namespace isaac
{
namespace bgzf
{

/**
 ** \brief Exception thrown when a zlib inflate method invocation fails.
 **
 **/
class BgzfInflateException: public common::IsaacException
{
public:
    BgzfInflateException(int error, z_stream &strm) : IsaacException(EINVAL, strm.msg ? strm.msg : "unknown error " + boost::lexical_cast<std::string>(error))
    {

    }
};


class BgzfReader
{
    // more blocks per pass reduces the amount of thread synchronization
    // unfortunately, the target buffer required gets too bit and
    // trashes the L3 cache. So, keep 1 for now
    static const unsigned BLOCKS_AT_ONCE = 1;
    common::FiniteCapacityVector<char, BLOCKS_AT_ONCE * 0xFFFF> compressedBlockBuffer_;

    z_stream strm_;
    char zalbuffer[3][65535];
    unsigned zalbufferFree;

public:
    BgzfReader() : zalbufferFree(sizeof(zalbuffer) / sizeof(zalbuffer[0]))
    {
        memset(&strm_, 0, sizeof(strm_));
    }
    unsigned readNextBlock(std::istream &is);
    void uncompressCurrentBlock(char* p, const std::size_t size);

private:
    void reset()
    {
        if (strm_.state)
        {
            inflateEnd(&strm_);
        }
        strm_.zalloc = zalloc;
        strm_.zfree = zfree;
        strm_.opaque = this;
        strm_.avail_in = 0;
        strm_.next_in = Z_NULL;
        int ret = inflateInit2(&strm_, 16+MAX_WBITS);
        if (Z_OK != ret)
        {
            BOOST_THROW_EXCEPTION(BgzfInflateException(ret, strm_));
        }
    }
    static voidpf zalloc OF((voidpf opaque, uInt items, uInt size))
    {
        ISAAC_ASSERT_MSG(reinterpret_cast<BgzfReader*>(opaque)->zalbufferFree, "Unexpected too many zalloc calls");
        ISAAC_ASSERT_MSG(sizeof(reinterpret_cast<BgzfReader*>(opaque)->zalbuffer[reinterpret_cast<BgzfReader*>(opaque)->zalbufferFree]) >= items * size, "Unexpected buffer size passed to zalloc : size=" << size << " items=" << items);
        --reinterpret_cast<BgzfReader*>(opaque)->zalbufferFree;
        return reinterpret_cast<BgzfReader*>(opaque)->zalbuffer[reinterpret_cast<BgzfReader*>(opaque)->zalbufferFree];
    }
    static void zfree OF((voidpf opaque, voidpf address))
    {
//        ISAAC_ASSERT_MSG(!reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbufferFree, "Unexpected unmatched zfree call");
//        ISAAC_ASSERT_MSG(reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbuffer == address, "Unexpected address passed to zfree");
        ++reinterpret_cast<BgzfReader*>(opaque)->zalbufferFree;
    }
};

class ParallelBgzfReader
{
    io::FileBufWithReopen fileBuffer_;
    std::istream is_;
    common::ThreadVector &threads_;
    const unsigned coresMax_;
    std::vector<BgzfReader> readers_;
    /// -1UL or the absolute offset of uncompressed data available in a threadBlocks_
    std::vector<unsigned long> threadOffsets_;
    /// uncompressed block of data which belongs at the corresponding threadOffsets_
    std::vector<unsigned long> threadBlockSizes_;
    /// size of the block that did not fit in the last result buffer
    unsigned long pendingBlockSize_;
    /// offset where the next bam block should decompress
    unsigned long nextUncompressedOffset_;

    boost::mutex stateMutex_;
    boost::condition_variable stateChangedCondition_;
    bool loadSlotAvailable_;
    unsigned computeSlotsAvailable_;

public:
    ParallelBgzfReader(
        common::ThreadVector &threads,
        const unsigned coresMax):
        fileBuffer_(std::ios_base::binary|std::ios_base::in),
        is_(&fileBuffer_),
        threads_(threads),
        coresMax_(coresMax),
        readers_(coresMax_),
        threadOffsets_(readers_.size(), -1UL),
        threadBlockSizes_(readers_.size(), -1UL),
        pendingBlockSize_(0),
        nextUncompressedOffset_(0),
        loadSlotAvailable_(true),
        computeSlotsAvailable_(readers_.size())
    {
    }
    void open(const boost::filesystem::path &bamPath);
    bool readMoreData(std::vector<char> &buffer);

    //    bool isOpen() const {return fileBuffer_.is_open();}
    bool isEof() {return is_.eof();}
private:
    void readMoreDataParallel(const unsigned threadNumber, std::vector<char> &buffer);

    void waitForLoadSlot(boost::unique_lock<boost::mutex> &lock);
    void releaseLoadSlot();
    void waitForComputeSlot(boost::unique_lock<boost::mutex> &lock);
    void releaseComputeSlot();
};

} // namespace bgzf
} // namespace isaac

#endif // #ifndef iSAAC_BGZF_BGZF_READER_HH
