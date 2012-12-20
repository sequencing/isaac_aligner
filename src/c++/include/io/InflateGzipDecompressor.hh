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
 ** \file InflateGzipDecompressor.hh
 **
 ** Decompression helper that goes directly to zlib inflate
 **
 ** \author Roman Petrovski
 **/

#include <zlib.h>

#ifndef iSAAC_IO_INFLATE_GZIP_DECOMPRESSOR_HH
#define iSAAC_IO_INFLATE_GZIP_DECOMPRESSOR_HH

//#pragma GCC push_options
//#pragma GCC optimize ("0")

namespace isaac
{
namespace io
{

/**
 ** \brief Exception thrown when a libxslt method invocation fails.
 **
 **/
class ZlibInflateException: public common::IsaacException
{
public:
    ZlibInflateException(int error, z_stream &strm) : IsaacException(EINVAL, strm.msg ? strm.msg : "unknown error " + boost::lexical_cast<std::string>(error))
    {

    }
};

/**
 * \brief Same as boost::gzip::decompressor, but ensures zlib buffer is allocated during construction.
 */
template <typename ContainterT>
class InflateGzipDecompressor
{
    ContainterT temporaryBuffer_;
    std::streamsize pendingBytes_;
    z_stream strm_;
    char zalbuffer[3][65535];
    unsigned zalbufferFree;
public:
    InflateGzipDecompressor() :
        pendingBytes_(0),
        zalbufferFree(sizeof(zalbuffer) / sizeof(zalbuffer[0]))
    {
        memset(&strm_, 0, sizeof(strm_));
    }

    InflateGzipDecompressor(const std::streamsize maxBufferSize) :
        pendingBytes_(0),
        zalbufferFree(sizeof(zalbuffer) / sizeof(zalbuffer[0]))
    {
        memset(&strm_, 0, sizeof(strm_));
        resize(maxBufferSize);
    }

    InflateGzipDecompressor(const InflateGzipDecompressor &that) :
        pendingBytes_(0),
        zalbufferFree(sizeof(zalbuffer) / sizeof(zalbuffer[0]))
    {
        memset(&strm_, 0, sizeof(strm_));
        resize(that.temporaryBuffer_.size());
    }

    void resize(const std::streamsize maxBufferSize)
    {
        temporaryBuffer_.resize(maxBufferSize);
    }

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
            BOOST_THROW_EXCEPTION(ZlibInflateException(ret, strm_));
        }

//        ISAAC_THREAD_CERR << "InflateGzipDecompressor::doSomethingToGetZlibInflateToAllocateWindowBuffer " << std::endl;
    }

    /**
     * \brief Same effect as gzip_decompresor::read but avoids dynamic memory allocation by
     *        using the externally supplied buffer.
     */
    std::streamsize read(
        std::istream &compressedStream,
        char *resultBuffer,
        std::streamsize resultBufferSize)
    {
        strm_.next_out = reinterpret_cast<Bytef *>(resultBuffer);
        strm_.avail_out = resultBufferSize;
        if (pendingBytes_)
        {
            strm_.next_in = reinterpret_cast<Bytef *>(&temporaryBuffer_.front());
            strm_.avail_in = pendingBytes_;
//            ISAAC_THREAD_CERR << "InflateGzipDecompressor::read writing " << pendingBytes_ <<
//                " pending bytes. " << std::endl;
            int err = inflate(&strm_, Z_SYNC_FLUSH);
//            const std::streamsize uncompressed = resultBufferSize - strm_.avail_out;
//            ISAAC_THREAD_CERR << "InflateGzipDecompressor::read inflate err: " << err << " uncompressed: " << uncompressed << std::endl;
            if (Z_OK != err && Z_STREAM_END != err)
            {
                BOOST_THROW_EXCEPTION(ZlibInflateException(err, strm_));
            }
            const std::streamsize compressedWritten = pendingBytes_ - strm_.avail_in;
            std::copy(temporaryBuffer_.begin() + compressedWritten,
                      temporaryBuffer_.begin() + pendingBytes_,
                      temporaryBuffer_.begin());
            pendingBytes_ = strm_.avail_in;

//            ISAAC_THREAD_CERR << "InflateGzipDecompressor::read written " << compressedWritten <<
//                " pending bytes. pendingBytes_: " << pendingBytes_ << std::endl;
        }
        if (strm_.avail_out)
        {
            if (compressedStream.good())
            {
                compressedStream.read(&temporaryBuffer_.front() + pendingBytes_, temporaryBuffer_.size() - pendingBytes_);
                pendingBytes_ += compressedStream.gcount();
            }
            if (!compressedStream.good() && !compressedStream.eof())
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to read compressed data"));
            }

            if (pendingBytes_)
            {
                strm_.next_in = reinterpret_cast<Bytef *>(&temporaryBuffer_.front());
                strm_.avail_in = pendingBytes_;
//                ISAAC_THREAD_CERR << "InflateGzipDecompressor::read writing " << pendingBytes_ <<
//                    " read bytes. " << std::endl;
                int err = inflate(&strm_, Z_SYNC_FLUSH);
//                const std::streamsize uncompressed = resultBufferSize - strm_.avail_out;
//                ISAAC_THREAD_CERR << "InflateGzipDecompressor::read inflate err: " << err << " uncompressed: " << uncompressed << std::endl;
                if (Z_OK != err && Z_STREAM_END != err)
                {
                    BOOST_THROW_EXCEPTION(ZlibInflateException(err, strm_));
                }
                const std::streamsize compressedWritten = pendingBytes_ - strm_.avail_in;
                std::copy(temporaryBuffer_.begin() + compressedWritten,
                          temporaryBuffer_.begin() + pendingBytes_,
                          temporaryBuffer_.begin());
                pendingBytes_ = strm_.avail_in;

//                ISAAC_THREAD_CERR << "InflateGzipDecompressor::read written " << compressedWritten <<
//                    " read bytes. pendingBytes_: " << pendingBytes_ << std::endl;
            }
        }
        std::streamsize ret = resultBufferSize - strm_.avail_out;

        if (!ret)
        {
            ISAAC_ASSERT_MSG(!compressedStream.good(), "When no bytes come out of decompressor expecting the input stream to be over");
            ISAAC_ASSERT_MSG(!pendingBytes_, "When no bytes come out of decompressor and "
                "the input stream is all finished, expecting the pendingBytes_ to be 0. Actual: " << pendingBytes_);

            ISAAC_THREAD_CERR << "InflateGzipDecompressor::read finished " << std::endl;

//            // NOTE: it is important that the resultBuffer can fit all the pending output.
//            close(sink, BOOST_IOS::out);
        }

        if (!ret)
        {
            // boost iostream filter must return -1 when no more data is available
            return -1;
        }

        return ret;
    }
private:
    static voidpf zalloc OF((voidpf opaque, uInt items, uInt size))
    {
        ISAAC_ASSERT_MSG(reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbufferFree, "Unexpected too many zalloc calls");
        ISAAC_ASSERT_MSG(sizeof(reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbuffer[reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbufferFree]) >= items * size, "Unexpected buffer size passed to zalloc : size=" << size << " items=" << items);
        --reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbufferFree;
        return reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbuffer[reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbufferFree];
    }
    static void zfree OF((voidpf opaque, voidpf address))
    {
//        ISAAC_ASSERT_MSG(!reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbufferFree, "Unexpected unmatched zfree call");
//        ISAAC_ASSERT_MSG(reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbuffer == address, "Unexpected address passed to zfree");
        ++reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbufferFree;
    }

};


} // namespace io
} // namespace isaac

//#pragma GCC pop_options

#endif // #ifndef iSAAC_IO_INFLATE_GZIP_DECOMPRESSOR_HH
