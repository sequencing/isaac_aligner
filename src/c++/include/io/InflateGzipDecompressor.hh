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
 * \brief serialization helper in this namespace to allow the boost internals to find the serialization operator
 *        without placing it in global namespace
 */
struct z_stream_serialization : public z_stream
{
    z_stream_serialization(const z_stream &that) : z_stream(that){}
};

inline std::ostream &operator <<(std::ostream &os, const z_stream_serialization &zs)
{
    return os << "z_stream(" <<
        " next_in:" << (void*)zs.next_in <<
        " avail_in:" << zs.avail_in <<
        " total_in:" << zs.total_in <<
        " next_out:" << (void*)zs.next_out <<
        " avail_out:" << zs.avail_out <<
        " total_out:" << zs.total_out <<
        " msg:" << (zs.msg ? (const char*)zs.msg : (const char*)("null"))<<
        " state:" << (void*)zs.state <<
        " zalloc:" << zs.zalloc <<
        " zfree:" << zs.zfree <<
        " opaque:" << zs.opaque <<
        " data_type:" << zs.data_type <<
        " adler:" << zs.adler <<
        " reserved:" << zs.reserved <<
        ")";
}

/**
 ** \brief Exception thrown when a zlib inflate method invocation fails.
 **
 **/
class ZlibInflateException: public common::IsaacException
{
public:
    ZlibInflateException(int error, z_stream &strm, const char *msg = 0) :
        IsaacException(EINVAL,
                       (strm.msg ?
                               strm.msg :
                               msg ?
                                   std::string(msg) + " unknown error " + boost::lexical_cast<std::string>(error) :
                                   "unknown error " + boost::lexical_cast<std::string>(error)) +
                                   boost::lexical_cast<std::string>(z_stream_serialization(strm)))

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
    static const std::size_t ALLOCATIONS_MAX = 3;
    static const std::size_t DEFAULT_BUFFER_SIZE = 4096 * 8;
    char zalbuffer[ALLOCATIONS_MAX][65535];
    unsigned zalbufferFree;
    bool allocationsBlocked_;
public:
    InflateGzipDecompressor() :
        pendingBytes_(0),
        zalbufferFree(ALLOCATIONS_MAX),
        allocationsBlocked_(false)
    {
        memset(&strm_, 0, sizeof(strm_));
        resize(DEFAULT_BUFFER_SIZE);
    }

    InflateGzipDecompressor(const std::streamsize maxBufferSize) :
        pendingBytes_(0),
        zalbufferFree(ALLOCATIONS_MAX),
        allocationsBlocked_(false)
    {
        memset(&strm_, 0, sizeof(strm_));
        resize(maxBufferSize);
    }

    InflateGzipDecompressor(const InflateGzipDecompressor &that) :
        pendingBytes_(0),
        zalbufferFree(ALLOCATIONS_MAX),
        allocationsBlocked_(false)
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
        clearStreamState();
        pendingBytes_ = 0;
    }

    /**
     * \brief Allows skipping some of the uncompressed data. Also supports concatenation of gzipped blocks
     *        by trying to read resultBufferSize bytes.
     */
    std::streamsize read(
        std::istream &compressedStream,
        std::streamsize skipUncompressedBytes,
        char *resultBuffer,
        std::streamsize resultBufferSize)
    {
        while (skipUncompressedBytes)
        {
            const std::streamsize uncompressed =
                read(compressedStream, resultBuffer, std::min(skipUncompressedBytes, resultBufferSize));
            if (-1 == uncompressed)
            {
                return -1;
            }
            ISAAC_ASSERT_MSG(uncompressed, "Uncompressed " << uncompressed << " bytes");
            skipUncompressedBytes -= uncompressed;
        }
        std::streamsize uncompressedTotal = 0;
        while(resultBufferSize)
        {
            const std::streamsize uncompressed = read(compressedStream, resultBuffer, resultBufferSize);
            if (-1 == uncompressed)
            {
                return uncompressedTotal ? uncompressedTotal : -1;
            }
            uncompressedTotal += uncompressed;
            resultBuffer += uncompressed;
            resultBufferSize -= uncompressed;
        }
        return uncompressedTotal;
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
        std::streamsize ret = 0;
        while (!ret)
        {
            if (pendingBytes_)
            {
                pendingBytes_ = processPendingBytes(strm_, pendingBytes_, temporaryBuffer_, !compressedStream.good());
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
                    pendingBytes_ = processPendingBytes(strm_, pendingBytes_, temporaryBuffer_, !compressedStream.good());
                }
            }

            ret = resultBufferSize - strm_.avail_out;
            if (!ret)
            {
                if (!pendingBytes_ && !compressedStream.good())
                {
                    break;
                }
                // assume concatenation of gz files. Reset zlib and start over
                resetStreamState();
            }
        }

        if (!ret)
        {
            // boost iostream filter must return -1 when no more data is available
            return -1;
        }

        return ret;
    }
private:
    /**
     * \brief Completely clears stream state. Resets input and output buffer information
     */
    void clearStreamState()
    {
        if (strm_.state)
        {
            inflateEnd(&strm_);
        }
        memset(&strm_, 0, sizeof(strm_));
        initStreamState();
    }
    /**
     * \brief Reset the stream state so that new compressed block can be processed. Does not touch input or output positions
     */
    void resetStreamState()
    {
        if (strm_.state)
        {
            inflateEnd(&strm_);
        }
        initStreamState();
    }

    /**
     * \brief Initializes decompression stream
     */
    void initStreamState()
    {
        strm_.zalloc = zalloc;
        strm_.zfree = zfree;
        strm_.opaque = this;

        int ret = inflateInit2(&strm_, 16+MAX_WBITS);
        if (Z_OK != ret)
        {
            BOOST_THROW_EXCEPTION(ZlibInflateException(ret, strm_));
        }
    }

    static std::streamsize processPendingBytes(z_stream &strm, const std::streamsize pendingBytes, ContainterT &temporaryBuffer, const bool endOfData)
    {
        strm.next_in = reinterpret_cast<Bytef *>(&temporaryBuffer.front());
        strm.avail_in = pendingBytes;

        int err = inflate(&strm, Z_SYNC_FLUSH);

        if (Z_OK != err && Z_STREAM_END != err)
        {
            if (endOfData)
            {
                BOOST_THROW_EXCEPTION(ZlibInflateException(err, strm, "Premature end of compressed stream reached. "));
            }
            else
            {
                BOOST_THROW_EXCEPTION(ZlibInflateException(err, strm));
            }
        }

        const std::streamsize compressedWritten = pendingBytes - strm.avail_in;
        std::copy(temporaryBuffer.begin() + compressedWritten,
                  temporaryBuffer.begin() + pendingBytes,
                  temporaryBuffer.begin());
        return strm.avail_in;
    }
    static voidpf zalloc OF((voidpf opaque, uInt items, uInt size))
    {
        ISAAC_ASSERT_MSG(!reinterpret_cast<InflateGzipDecompressor*>(opaque)->allocationsBlocked_, "TODO: implement sensible memory management here");
        ISAAC_ASSERT_MSG(reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbufferFree, "Unexpected too many zalloc calls");
        ISAAC_ASSERT_MSG(sizeof(reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbuffer[reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbufferFree]) >= items * size, "Unexpected buffer size passed to zalloc : size=" << size << " items=" << items);
        --reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbufferFree;
        return reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbuffer[reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbufferFree];
    }
    static void zfree OF((voidpf opaque, voidpf address))
    {
        // this is the memory management for poor. The assumption is that there will be no more than ALLOCATIONS_MAX zalloc calls and that
        // ALLOCATIONS_MAX zfree calls will occur without any zalloc calls in between them.
        reinterpret_cast<InflateGzipDecompressor*>(opaque)->allocationsBlocked_ =
            ALLOCATIONS_MAX != ++reinterpret_cast<InflateGzipDecompressor*>(opaque)->zalbufferFree;
    }

};


} // namespace io
} // namespace isaac

//#pragma GCC pop_options

#endif // #ifndef iSAAC_IO_INFLATE_GZIP_DECOMPRESSOR_HH
