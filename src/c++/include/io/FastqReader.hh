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
 ** \file FastqReader.hh
 **
 ** Component to read FASTQ files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_FASTQ_READER_HH
#define iSAAC_IO_FASTQ_READER_HH

#include <vector>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/FiniteCapacityVector.hh"
#include "flowcell/ReadMetadata.hh"
#include "io/InflateGzipDecompressor.hh"
#include "io/FileBufCache.hh"
#include "oligo/Nucleotides.hh"

namespace isaac
{
namespace io
{

/**
 ** \brief Exception thrown when fastq violation is encountered.
 **
 **/
class FastqFormatException: public common::IsaacException
{
public:
    FastqFormatException(const std::string &message):
        common::IsaacException(EINVAL, message){}
};

class FastqReader: boost::noncopyable
{
    static const std::size_t UNCOMPRESSED_BUFFER_SIZE = 4096;
    // Since we're inflating, it is important that the uncompressed buffer does not
    // run out of space for each individual chunk of compressed bytes pushed into decompressor.
    static const std::size_t DECOMPRESSOR_BUFFER_SIZE = UNCOMPRESSED_BUFFER_SIZE;
    static const unsigned FASTQ_QSCORE_OFFSET = 33;
    const bool allowVariableLength_;

    FileBufWithReopen fileBuffer_;
    typedef common::FiniteCapacityVector<char, UNCOMPRESSED_BUFFER_SIZE> BufferType;
    io::InflateGzipDecompressor<BufferType> decompressor_;

    //boost::filesystem::path forces intermediate string construction during reassignment...
    std::string fastqPath_;
    bool compressed_;
    bool reachedEof_;
    std::size_t filePos_;

    BufferType buffer_;
    BufferType::const_iterator headerBegin_;
    BufferType::const_iterator headerEnd_;
    BufferType::const_iterator baseCallsBegin_;
    BufferType::const_iterator baseCallsEnd_;
    BufferType::const_iterator qScoresBegin_;
    BufferType::const_iterator endIt_;

    static const oligo::Translator translator_;

public:
    static const unsigned INCORRECT_FASTQ_BASE = 5;

    FastqReader(const bool allowVariableLength);
    FastqReader(const bool allowVariableLength, const boost::filesystem::path &fastqPath);

    void reservePathBuffers(std::size_t maxPathLength)
    {
        fastqPath_ = std::string(maxPathLength, ' ').c_str();
    }

    void open(const boost::filesystem::path &fastqPath);

    void next();

    template <typename InsertIt>
    void getBcl(const flowcell::ReadMetadata &readMetadata, InsertIt &it) const;

    const std::string &getPath() const
    {
        return fastqPath_;
    }

    std::size_t getRecordOffset() const
    {
        return getOffset(headerBegin_);
    }

    bool hasData() const
    {
        return (!reachedEof_ || buffer_.end() != headerBegin_);
    }

    typedef std::pair<BufferType::const_iterator, BufferType::const_iterator > IteratorPair;
    IteratorPair getHeader() const
    {
        return std::make_pair(headerBegin_, headerEnd_);
    }

    unsigned getReadLength() const
    {
        return std::distance(baseCallsBegin_, baseCallsEnd_);
    }

private:
    typedef boost::error_info<struct tag_errmsg, std::string> errmsg_info;

    void resetBuffer();
    std::size_t getOffset(BufferType::const_iterator it) const;
    void findHeader();
    void findSequence();
    void findQScores();
    void findQScoresEnd();
    bool fetchMore();

    std::size_t readCompressedFastq(std::istream &is, char *buffer, std::size_t amount);
    std::size_t readFlatFastq(std::istream &is, char *buffer, std::size_t amount);
};

template <typename InsertIt>
void FastqReader::getBcl(const flowcell::ReadMetadata &readMetadata, InsertIt &it) const
{
    const InsertIt start = it;
    BufferType::const_iterator baseCallsIt = baseCallsBegin_;
    BufferType::const_iterator qScoresIt = qScoresBegin_;
    std::vector<unsigned>::const_iterator cycleIterator = readMetadata.getCycles().begin();
    unsigned currentCycle = readMetadata.getFirstReadCycle();
    for(;endIt_ != qScoresIt && readMetadata.getCycles().end() != cycleIterator; ++baseCallsIt, ++qScoresIt, ++currentCycle)
    {
//        ISAAC_THREAD_CERR << "cycle " << *cycleIterator << std::endl;
        if (*cycleIterator != currentCycle)
        {
            continue;
        }
        const unsigned char baseValue = translator_[*baseCallsIt];
        if (oligo::invalidOligo == baseValue)
        {
            *it = 0;
        }
        else if(INCORRECT_FASTQ_BASE == baseValue)
        {
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format("Invalid oligo %c found in %s at offset %u") %
                *baseCallsIt % getPath() % getOffset(baseCallsIt)).str()));
        }
        else
        {
            const unsigned char baseQuality = (*qScoresIt - FASTQ_QSCORE_OFFSET);
            if ((1 << 6) <= baseQuality)
            {
                BOOST_THROW_EXCEPTION(FastqFormatException((boost::format("Invalid quality %d found in %s at offset %u. Base quality scores [0-63] supported only.") %
                    unsigned(baseQuality) % getPath() % getOffset(baseCallsIt)).str()));
            }
            *it = baseValue | (baseQuality << 2);
        }
//        ISAAC_THREAD_CERR << "stored cycle " << *cycleIterator << std::endl;
        ++it;
        ++cycleIterator;
    }

    const unsigned extracted = std::distance<std::vector<unsigned>::const_iterator>(readMetadata.getCycles().begin(), cycleIterator);
    if (!allowVariableLength_)
    {
        if (readMetadata.getCycles().size() != extracted)
        {
    //            ISAAC_ASSERT_MSG(false, "Fastq read shorter than expected." << ret);
            BOOST_THROW_EXCEPTION(common::IoException(EINVAL, (boost::format("Read length (%d) "
                " is different from expected %d in %s:%u. Record %s") %
                extracted % readMetadata.getCycles().size() %
                getPath() % getRecordOffset() %
                std::string(headerBegin_, endIt_)).str()));
        }
        else
        {
//            ISAAC_THREAD_CERR << "length ok " << std::endl;30=.cilprw
        }
    }
    else if (extracted < readMetadata.getCycles().size())
    {
        it = std::fill_n(it, readMetadata.getCycles().size() - extracted, 0);
    }

    ISAAC_ASSERT_MSG(readMetadata.getCycles().size() == std::size_t(std::distance(start, it)),
        "unexpected number of cycles read: " << std::distance(start, it) << " expected: " << readMetadata);
}
} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_FASTQ_READER_HH
