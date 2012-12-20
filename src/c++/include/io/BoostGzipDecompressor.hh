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
 ** \file BoostGzipDecompressor.hh
 **
 ** Same as boost::gzip::decompressor, but ensures zlib buffer is allocated during construction.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_BOOST_GZIP_DECOMPRESSOR_HH
#define iSAAC_IO_BOOST_GZIP_DECOMPRESSOR_HH

#include <boost/format.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace isaac
{
namespace io
{

/**
 * \brief Tiny helper class for appending decompressor output to a fixed size buffer
 */
class BufferInsertDevice {
public:
    typedef char  char_type;
    typedef boost::iostreams::sink_tag   category;
    BufferInsertDevice(char *pData, unsigned size) : pData_(pData), sizeLeft_(size) { }
    std::streamsize write(const char_type* s, std::streamsize n)
    {
        //ISAAC_ASSERT_MSG(sizeLeft_ >= n, "Buffer overflow");
        std::streamsize toCopy = std::min(n, sizeLeft_);
        pData_ = std::copy(s, s + toCopy, pData_);
        sizeLeft_ -= toCopy;
        return toCopy;
    }

    std::streamsize getSizeLeft() const {return sizeLeft_;}
protected:
    char *pData_;
    std::streamsize sizeLeft_;
};

/**
 * \brief Same as boost::gzip::decompressor, but ensures zlib buffer is allocated during construction.
 */
template <typename ContainterT>
class BoostGzipDecompressor : public boost::iostreams::gzip_decompressor
{
    ContainterT temporaryBuffer_;
    std::streamsize pendingBytes_;
public:
    BoostGzipDecompressor() :
        boost::iostreams::gzip_decompressor(),
        pendingBytes_(0)
    {
        doSomethingToGetZlibInflateToAllocateWindowBuffer();
    }

    BoostGzipDecompressor(const std::streamsize maxBufferSize) :
        boost::iostreams::gzip_decompressor(),
        pendingBytes_(0)
    {
        temporaryBuffer_.resize(maxBufferSize);
        doSomethingToGetZlibInflateToAllocateWindowBuffer();
    }

    BoostGzipDecompressor(const BoostGzipDecompressor &that) :
        boost::iostreams::gzip_decompressor(that),
        pendingBytes_(0)
    {
        temporaryBuffer_.resize(that.temporaryBuffer_.capacity());
        doSomethingToGetZlibInflateToAllocateWindowBuffer();
    }

    void reserve(const std::streamsize maxBufferSize)
    {
        temporaryBuffer_.resize(maxBufferSize);
    }

    /**
     * \brief This Voodoo serves only one purpose: get to the point where zlib inflate is
     *        called so that zlib allocates window buffer during the decompressor
     *        construction and not at the runtime. TODO: If you known a cleaner way to do it, please improve.
     */
    void doSomethingToGetZlibInflateToAllocateWindowBuffer()
    {
        boost::iostreams::filtering_istream stm;

        static const char uncompressedStuff[] = "t";
        boost::iostreams::array_source uncompSrc(uncompressedStuff, sizeof(uncompressedStuff));

        // there will be some room needed for header
        char compressedStuff[1024];
        boost::iostreams::gzip_compressor comp;
        stm.push(uncompSrc);
        unsigned compressed = comp.read(stm, compressedStuff, sizeof(compressedStuff));

        stm.pop();
        boost::iostreams::array_source compSrc(compressedStuff, compressed);

        stm.push(compSrc);
        char uncompressedStuff2[sizeof(uncompressedStuff)];
        boost::iostreams::gzip_decompressor::read(stm, uncompressedStuff2, 1);
        close(stm, BOOST_IOS::in);
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
        io::BufferInsertDevice sink(resultBuffer, resultBufferSize);
        if (pendingBytes_)
        {
//            ISAAC_THREAD_CERR << "BoostGzipDecompressor::read writing " << pendingBytes_ <<
//                " pending bytes. " << std::endl;
            const std::streamsize compressedWritten = write(sink, &temporaryBuffer_.front(), pendingBytes_);
            std::copy(temporaryBuffer_.begin() + compressedWritten,
                      temporaryBuffer_.begin() + pendingBytes_,
                      temporaryBuffer_.begin());
            pendingBytes_ -= compressedWritten;

//            ISAAC_THREAD_CERR << "BoostGzipDecompressor::read written " << compressedWritten <<
//                " pending bytes. pendingBytes_: " << pendingBytes_ << std::endl;
        }
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
//            ISAAC_THREAD_CERR << "BoostGzipDecompressor::read writing " << pendingBytes_ <<
//                " read bytes. " << std::endl;
            const std::streamsize compressedWritten = write(sink, &temporaryBuffer_.front(), pendingBytes_);
            std::copy(temporaryBuffer_.begin() + compressedWritten,
                      temporaryBuffer_.begin() + pendingBytes_,
                      temporaryBuffer_.begin());
            pendingBytes_ -= compressedWritten;

//            ISAAC_THREAD_CERR << "BoostGzipDecompressor::read written " << compressedWritten <<
//                " read bytes. pendingBytes_: " << pendingBytes_ << std::endl;
        }
        std::streamsize ret = resultBufferSize - sink.getSizeLeft();
        if (!ret)
        {
            ISAAC_ASSERT_MSG(!compressedStream.good(), "When no bytes come out of decompressor expecting the input stream to be over");
            ISAAC_ASSERT_MSG(!pendingBytes_, "When no bytes come out of decompressor and "
                "the input stream is all finished, expecting the pendingBytes_ to be 0. Actual: " << pendingBytes_);

            ISAAC_THREAD_CERR << "BoostGzipDecompressor::read finished " << std::endl;
        }

        if (!ret)
        {
            // boost iostream filter must return -1 when no more data is available
            return -1;
        }

        return ret;
    }
};


} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_BOOST_GZIP_DECOMPRESSOR_HH
