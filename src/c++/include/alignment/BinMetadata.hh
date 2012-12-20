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
 ** \file BinMetadata.hh
 **
 ** \brief Metadata associated to the unsorted alignment results
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_BIN_METADATA_HH
#define iSAAC_ALIGNMENT_BIN_METADATA_HH

#include <numeric>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "flowcell/BarcodeMetadata.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace alignment
{

struct BinMetadata;
inline std::ostream &operator<<(std::ostream &os, const BinMetadata &binMetadata);

struct BarcodeCounts
{
    BarcodeCounts() : elements_(0), gaps_(0), cigarLength_(0){}
    // total number of elements in the bin barcode
    unsigned long elements_;
    // total number of gaps in the bin barcode
    unsigned long gaps_;
    // sum of all fragment cigar lengths in the bin barcode.
    unsigned long cigarLength_;
};


struct BinChunk
{
    std::vector<BarcodeCounts> barcodeBreakdown_;
    unsigned long dataSize_;

    // required for serialization
    BinChunk(): dataSize_(0){}
    BinChunk(const unsigned barcodesCount) : barcodeBreakdown_(barcodesCount), dataSize_(0){}

    unsigned long getTotalCigarLength() const
    {
        return std::accumulate(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), 0,
                               boost::bind(std::plus<unsigned long>(),
                                           _1, boost::bind(&BarcodeCounts::cigarLength_, _2)));
    }

    unsigned long getTotalElements() const
    {
        return std::accumulate(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), 0,
                               boost::bind(std::plus<unsigned long>(),
                                           _1, boost::bind(&BarcodeCounts::elements_, _2)));
    }

    unsigned long getBarcodeGapCount(const unsigned barcodeIdx) const
        {return barcodeBreakdown_.at(barcodeIdx).gaps_;}

    unsigned long getBarcodeElements(const unsigned barcodeIdx) const
        {return barcodeBreakdown_.at(barcodeIdx).elements_;}

    void incrementGapCount(const unsigned long by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).gaps_ += by;}

    void incrementCigarLength(const unsigned long by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).cigarLength_ += by;}

    void incrementSeIdxElements(const unsigned long by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).elements_ += by;}

    void incrementRIdxElements(const unsigned long by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).elements_ += by;}

    void incrementFIdxElements(const unsigned long by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).elements_ += by;}

    void incrementNmElements(const unsigned long by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).elements_ += by;}
};

class BinDataDistribution : std::vector<BinChunk>
{
    unsigned long chunkSize_;
    bool offsetsTallied_;

public:
    using std::vector<BinChunk>::reserve;
    using std::vector<BinChunk>::size;
    BinDataDistribution(const unsigned barcodesCount, unsigned long length)
    // one more chunk so that tallyOffset produces the end offset for the last present chunk
        :std::vector<BinChunk>(length / getChunkSize(length) + 2UL, BinChunk(barcodesCount)),
         chunkSize_(getChunkSize(length)), offsetsTallied_(false)
    {
    }

    unsigned long addBytes(unsigned long binGenomicOffset, unsigned count);
    unsigned long tallyOffsets();

    BinDataDistribution &operator =(const BinDataDistribution &that);

    BinChunk &getChunk(unsigned long binGenomicOffset)
        {return this->operator[](getIndex(binGenomicOffset));}

    std::size_t getIndex(unsigned long binGenomicOffset) const;
    unsigned long getChunkSize() const {return chunkSize_;}

    /// Aim to have 1024 chunks. This will require about 4096*1024 (4 megabytes) of cache when pre-sorting bin during the loading in bam generator
    static unsigned long getChunkSize(unsigned long length) {return length < 1024UL ? 1024UL : length / 1024UL;}

    unsigned long getChunkEndOffset(std::size_t chunk);

    unsigned long getTotalCigarLength() const;
    unsigned long getBarcodeGapCount(const unsigned barcodeIdx) const;
    unsigned long getBarcodeElements(const unsigned barcodeIdx) const;
    unsigned long getTotalElements() const;

    void incrementCigarLength(const unsigned long binGenomicOffset, const unsigned long by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementCigarLength(by, barcodeIdx);}

    void incrementGapCount(const unsigned long binGenomicOffset, const unsigned long by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementGapCount(by, barcodeIdx);}

    void incrementSeIdxElements(const unsigned long binGenomicOffset, const unsigned long by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementSeIdxElements(by, barcodeIdx);}

    void incrementRIdxElements(const unsigned long binGenomicOffset, const unsigned long by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementRIdxElements(by, barcodeIdx);}

    void incrementFIdxElements(const unsigned long binGenomicOffset, const unsigned long by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementFIdxElements(by, barcodeIdx);}

    void incrementNmElements(const unsigned long binGenomicOffset, const unsigned long by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementNmElements(by, barcodeIdx);}

    unsigned long removeChunksBefore(const unsigned long minOffset);
    unsigned long removeChunksAfter(const unsigned long minOffset);

    /*
     * \brief enable serialization
     */
    template <class Archive> friend void serialize(Archive &ar, BinDataDistribution &bm, const unsigned int version);
};

class BinMetadata
{
    unsigned binIndex_;
    /// first genomic position covered by the bin
    reference::ReferencePosition binStart_;
    /// bin length in bases
    unsigned long length_;
    boost::filesystem::path binFilePath_;
    boost::filesystem::path fIdxFilePath_;
    boost::filesystem::path rIdxFilePath_;
    boost::filesystem::path seIdxFilePath_;
    // offset from the beginning of the data file.
    // Note that single file can later be broken down into multiple BinMetadata objects
    unsigned long dataOffset_;
    // number of bytes stored in binFilePath_ at dataOffset_
    unsigned long dataSize_;
    unsigned long seIdxElements_;
    unsigned long rIdxElements_;
    unsigned long fIdxElements_;
    unsigned long nmElements_;

    BinDataDistribution dataDistribution_;

    /*
     * \brief enable serialization
     */
    template <class Archive> friend void serialize(Archive &ar, BinMetadata &bm, const unsigned int version);

public:
    BinMetadata() :
        binIndex_(0),
        binStart_(0),
        length_(0),
        dataOffset_(0),
        dataSize_(0),
        seIdxElements_(0),
        rIdxElements_(0),
        fIdxElements_(0),
        nmElements_(0),
        dataDistribution_(0,0){}

    BinMetadata(
        const unsigned barcodesCount,
        const unsigned binIndex,
        const reference::ReferencePosition binStart,
        const unsigned long length,
        const boost::filesystem::path &binFilepath) :
            binIndex_(binIndex),
            binStart_(binStart),
            length_(length),
            binFilePath_(binFilepath),
            fIdxFilePath_(getFIdxFilePath(binFilePath_)),
            rIdxFilePath_(getRIdxFilePath(binFilePath_)),
            seIdxFilePath_(getSeIdxFilePath(binFilePath_)),
            dataOffset_(0),
            dataSize_(0),
            seIdxElements_(0),
            rIdxElements_(0),
            fIdxElements_(0),
            nmElements_(0),
            dataDistribution_(barcodesCount, length_){}

    /**
     * \return BinMedata which guarantees to have the chunks with
     *         offset: minOffset <= offset < (minOffset + minSize)
     */
    BinMetadata getChunks(
        const unsigned long minOffset,
        const unsigned long minSize) const
    {
        ISAAC_ASSERT_MSG(isUnalignedBin(), "Splitting bins is supported only for unaligned bin");
        ISAAC_ASSERT_MSG(!rIdxElements_, "Splitting bins is supported only for unaligned bin");
        ISAAC_ASSERT_MSG(!fIdxElements_, "Splitting bins is supported only for unaligned bin");
        ISAAC_ASSERT_MSG(!seIdxElements_, "Splitting bins is supported only for unaligned bin");

        BinMetadata ret(*this);
        ret.removeChunksBefore(minOffset);
        ret.removeChunksAfter(minSize);

        return ret;
    }

    void removeChunksBefore(const unsigned long minOffset)
    {
        const unsigned long removedBytes = dataDistribution_.removeChunksBefore(minOffset);
        dataOffset_ += removedBytes;
        dataSize_ -= removedBytes;
    }

    void removeChunksAfter(const unsigned long minOffset)
    {
        dataSize_ = dataDistribution_.removeChunksAfter(minOffset);
    }

    unsigned getIndex() const
    {
        return binIndex_;
    }

    void setIndex(unsigned binIndex)
    {
        binIndex_ = binIndex;
    }

    reference::ReferencePosition getBinStart() const
    {
        return binStart_;
    }

    reference::ReferencePosition getBinEnd() const
    {
        return isUnalignedBin() ? reference::ReferencePosition(reference::ReferencePosition::NoMatch) : binStart_ + length_;
    }

    bool isUnalignedBin() const {return binStart_.isTooManyMatch();}

    const boost::filesystem::path & getPath() const
    {
        return binFilePath_;
    }
    const std::string &getPathString() const
    {
        return binFilePath_.string();
    }

    const boost::filesystem::path &getFIdxFilePath() const
    {
        return fIdxFilePath_;
    }

    const boost::filesystem::path &getRIdxFilePath() const
    {
        return rIdxFilePath_;
    }

    const boost::filesystem::path &getSeIdxFilePath() const
    {
        return seIdxFilePath_;
    }

    static boost::filesystem::path getFIdxFilePath(const boost::filesystem::path &binFilePath)
    {
        return binFilePath.string() + ".fw-idx";
    }

    static boost::filesystem::path getRIdxFilePath(const boost::filesystem::path &binFilePath)
    {
        return binFilePath.string() + ".rs-idx";
    }

    static boost::filesystem::path getSeIdxFilePath(const boost::filesystem::path &binFilePath)
    {
        return binFilePath.string() + ".se-idx";
    }

    unsigned long getDataOffset() const
    {
        return dataOffset_;
    }

    unsigned long getDataSize() const
    {
        return dataSize_;
    }

    unsigned long getLength() const {return length_;}

    void tallyChunkOffsets()
    {
        dataDistribution_.tallyOffsets();
        dataSize_ = 0;
    }

    /**
     * \brief increment the the corresponding chunk size and total data size.
     *
     * \return pair of total data size before increment and data offset of the corresponding chunk prior to being incremented
     */
    std::pair<unsigned long, unsigned long> incrementDataSize(const reference::ReferencePosition pos, const unsigned long by)
    {
        const std::pair<unsigned long, unsigned long> ret = std::make_pair(dataSize_, dataDistribution_.addBytes(getDataDistributionKey(pos), by));
        dataSize_ += by;
        return ret;
    }

    /**
     * \brief increment the the corresponding chunk size and total data size.
     *
     * \return pair of total data size before increment and data offset of the corresponding chunk prior to being incremented
     */
    std::pair<unsigned long, unsigned long> incrementDataSize(const unsigned long clusterNumber, const unsigned long by)
    {
        const std::pair<unsigned long, unsigned long> ret = std::make_pair(dataSize_, dataDistribution_.addBytes(getDataDistributionKey(clusterNumber), by));
        dataSize_ += by;
        return ret;
    }

    unsigned long getDataDistributionKey(const unsigned long clusterNumber)
    {
        ISAAC_ASSERT_MSG(isUnalignedBin(), "Aligned bins must use ReferencePosition as hash key." << *this);
        return clusterNumber;
    }

    unsigned long getDataDistributionKey(const reference::ReferencePosition pos)
    {
        ISAAC_ASSERT_MSG(!isUnalignedBin(), "Unaligned bins must use cluster number as key." << *this);
        ISAAC_ASSERT_MSG(pos >= binStart_, "pos: " << pos << " binStart_:" << binStart_);// << *this);
        ISAAC_ASSERT_MSG(pos < (binStart_ + length_), "pos: " << pos << " binStart_ + length_:" << (binStart_ + length_));//  << *this);
        return pos - binStart_;
    }

    unsigned long getSeIdxElements() const
    {
        return seIdxElements_;
    }

    void incrementSeIdxElements(const reference::ReferencePosition pos, const unsigned long by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementSeIdxElements(getDataDistributionKey(pos), by, barcodeIdx);
        seIdxElements_ += by;
    }

    unsigned long getRIdxElements() const
    {
        return rIdxElements_;
    }

    void incrementRIdxElements(const reference::ReferencePosition pos, const unsigned long by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementRIdxElements(getDataDistributionKey(pos), by, barcodeIdx);
        rIdxElements_ += by;
    }

    unsigned long getFIdxElements() const
    {
        return fIdxElements_;
    }

    void incrementFIdxElements(const reference::ReferencePosition pos, const unsigned long by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementFIdxElements(getDataDistributionKey(pos), by, barcodeIdx);
        fIdxElements_ += by;
    }

    unsigned long getNmElements() const
    {
        return nmElements_;
    }

    void incrementNmElements(const unsigned long sequenceHash, const unsigned long by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementNmElements(getDataDistributionKey(sequenceHash), by, barcodeIdx);
        nmElements_ += by;
    }

    void incrementGapCount(const reference::ReferencePosition pos, const unsigned long by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementGapCount(getDataDistributionKey(pos), by, barcodeIdx);
    }

    void incrementCigarLength(const reference::ReferencePosition pos, const unsigned long by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementCigarLength(getDataDistributionKey(pos), by, barcodeIdx);
    }

    unsigned long getTotalElements() const
    {
        //return  getSeIdxElements() + getRIdxElements() + getFIdxElements() + getNmElements();
        return dataDistribution_.getTotalElements();
    }

    unsigned long getBarcodeElements(const unsigned barcodeIdx) const
    {
        return dataDistribution_.getBarcodeElements(barcodeIdx);
    }

    unsigned long getBarcodeGapCount(const unsigned barcodeIdx) const
    {
        return dataDistribution_.getBarcodeGapCount(barcodeIdx);
    }

    unsigned long getTotalCigarLength() const
    {
        return dataDistribution_.getTotalCigarLength();
    }

    const BinDataDistribution &getDataDistribution() const
    {
        return dataDistribution_;
    }
};


struct BinMetadataList : public std::vector<alignment::BinMetadata>
{
    BinMetadataList(){}
    BinMetadataList(size_t size) : std::vector<alignment::BinMetadata>(size){}
};

inline std::ostream &operator<<(std::ostream &os, const BinMetadata &binMetadata)
{
    return os << "BinMetadata("
              << binMetadata.getIndex() << "id "
              << binMetadata.getBinStart() << "bs "
              << binMetadata.getLength() << "bl "
              << binMetadata.getDataSize() << "ds "
              << binMetadata.getDataOffset() << "do "
              << binMetadata.getSeIdxElements() << "se "
              << binMetadata.getRIdxElements() << "rs "
              << binMetadata.getFIdxElements() << "f "
              << binMetadata.getPathString() << ")";
}


/**
 * \param binOffset Offset relative to the beginning of the Bin
 *
 * \return data offset prior to incrementing
 */
inline unsigned long BinDataDistribution::addBytes(unsigned long binGenomicOffset, unsigned count)
{
    BinChunk &ref = getChunk(binGenomicOffset);
    const unsigned long ret = ref.dataSize_;
//        ISAAC_THREAD_CERR << "addBytes binOffset:" << binOffset << "getChunkSize()" << getChunkSize() <<
//            " count:" << count <<
//            " binIndex: " << binIndex << " ret:" << ret << std::endl;
    ref.dataSize_ += count;
    return ret;
}

/**
 * \brief replace contents of each chunk with sum of contents of all previous chunks
 *
 * \return the total number of bytes occupied by data
 */
inline unsigned long BinDataDistribution::tallyOffsets()
{
    unsigned long offset = 0;
    std::vector<BinChunk> &v = *this;
    BOOST_FOREACH(BinChunk &chunk, v)
    {
        using std::swap;
        swap(offset, chunk.dataSize_);
//            ISAAC_THREAD_CERR << "tallyOffsets: " << chunkOffset << std::endl;
        offset += chunk.dataSize_;
    }
    offsetsTallied_ = true;
    return offset;
}

inline BinDataDistribution &BinDataDistribution::operator =(const BinDataDistribution &that)
{
    assign(that.begin(), that.end());
    offsetsTallied_ = that.offsetsTallied_;
    return *this;
}

inline std::size_t BinDataDistribution::getIndex(unsigned long binGenomicOffset) const
{
    std::size_t ret = binGenomicOffset / getChunkSize();
    ISAAC_ASSERT_MSG(size() > ret, "bin index " << ret << " for offset " << binGenomicOffset << "is too big");
    return ret;
}

inline unsigned long BinDataDistribution::getChunkEndOffset(std::size_t chunk)
{
    ISAAC_ASSERT_MSG(offsetsTallied_, "getChunkEndOffset for untallied distribution");
    return size() > (chunk + 1) ? at(chunk+1).dataSize_ : back().dataSize_;
}

inline unsigned long BinDataDistribution::getTotalCigarLength() const
{
    return std::accumulate(begin(), end(), 0,
                           boost::bind(std::plus<unsigned long>(),
                                       _1, boost::bind(&BinChunk::getTotalCigarLength, _2)));
}

inline unsigned long BinDataDistribution::getBarcodeGapCount(const unsigned barcodeIdx) const
{
    return std::accumulate(begin(), end(), 0,
                           boost::bind(std::plus<unsigned long>(),
                                       _1, boost::bind(&BinChunk::getBarcodeGapCount, _2, barcodeIdx)));
}

inline unsigned long BinDataDistribution::getBarcodeElements(const unsigned barcodeIdx) const
{
    return std::accumulate(begin(), end(), 0,
                           boost::bind(std::plus<unsigned long>(),
                                       _1, boost::bind(&BinChunk::getBarcodeElements, _2, barcodeIdx)));
}

inline unsigned long BinDataDistribution::getTotalElements() const
{
    return std::accumulate(begin(), end(), 0,
                           boost::bind(std::plus<unsigned long>(),
                                       _1, boost::bind(&BinChunk::getTotalElements, _2)));
}

/**
 * \return number of bytes removed
 */
inline unsigned long BinDataDistribution::removeChunksBefore(const unsigned long minOffset)
{
    unsigned long offset = 0;
    std::vector<BinChunk>::iterator e = begin();
    for (; end() != e && minOffset > offset; ++e)
    {
        offset += e->dataSize_;
    }

    erase(begin(), e);
    return offset;
}

/**
 * \return number of bytes left
 */
inline unsigned long BinDataDistribution::removeChunksAfter(const unsigned long minOffset)
{
    unsigned long offset = 0;
    std::vector<BinChunk>::iterator b = begin();
    for (; end() != b && minOffset > offset; ++b)
    {
        offset += b->dataSize_;
    }

    erase(b, end());
    return offset;
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_BIN_METADATA_HH
