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
 ** \file BamIndexer.hh
 **
 ** \brief implements a boost iostreams filter that indexes a BAM input stream:
 ** forwards BAM stream to first output and generates a BAI index stream as second output
 **
 ** \author Lilian Janin
 **/

#ifndef iSAAC_BAM_BAM_INDEXER_HH
#define iSAAC_BAM_BAM_INDEXER_HH

#include <fstream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "common/Debug.hh"
#include "build/FragmentAccessorBamAdapter.hh"


namespace isaac
{
namespace bam
{

namespace bios=boost::iostreams;


// 512 Mbases is the longest chromosome length allowed in a BAM index
static const uint32_t BAM_MAX_CONTIG_LENGTH     = 512*1024*1024; 

// =(8^6-1)/7+1, as defined in samtools
static const uint32_t BAM_MAX_BIN               = 37450;         

// Each non-leaf bin contains 8 sub-bins => we expect a maximum of 7 clusters per bin, but we may sometimes
// get unlucky and a cluster may be split in two if some reads alternate between 2 bins just when they also
// reach the end of a BGZF block
static const uint32_t MAX_CLUSTER_PER_INDEX_BIN = 16;            

// BAM format constant
static const uint32_t BAM_FUNMAP = 4; 


class VirtualOffset
{
    uint64_t val_;

public:
    VirtualOffset() : val_(0) {}
    void set( uint64_t cOffset, uint32_t uOffset) { val_ = (cOffset << 16) | uOffset; }
    void set( uint64_t val)                           { val_ = val; }
    uint64_t get()                const               { return val_; }
    uint64_t compressedOffset()   const               { return val_>>16; }
    uint32_t uncompressedOffset() const               { return val_ & 0xFFFF; }

    friend std::ostream& operator<<( std::ostream& os, const VirtualOffset& virtualOffset );
};

inline std::ostream& operator<<( std::ostream& os, const VirtualOffset& virtualOffset )
{
    return os << "{" << (virtualOffset.val_ >> 16) << ", " << (virtualOffset.val_ & 0xFFFF) << "}";
}

typedef std::pair< VirtualOffset, VirtualOffset > VirtualOffsetPair;
typedef uint64_t UnresolvedOffset;


template<typename Device>
class BamIndexer
{
    // Single uncompressed BGZF chunks cannot contain more than 65535 bytes. Our uncompressed buffer contains
    // 1 uncompressed BGZF chunk plus the remainder of the previous BGZF chunk = 2 chunks in the worst case
    static const uint32_t MAX_UNCOMPRESSED_SIZE     = 65536*2;       

    // Single compressed BGZF chunk plus the remainder of the previous BGZF chunk = 2 chunks in the worst case.
    // Each compressed chunk is believed to be no larger than its uncompressed data, which is limited to 64KB
    static const uint32_t MAX_COMPRESSED_SIZE       = 65536*2;       

    // Buffer size for boost::gzip_decompressor, including some overhead for its internals
    static const uint32_t GZIP_INTERNAL_RAM = MAX_UNCOMPRESSED_SIZE * 2;

public:
    typedef char char_type;
    struct category : bios::multichar_output_filter_tag , bios::flushable_tag {};

    BamIndexer(Device& baiSink);
    BamIndexer(const BamIndexer& that);
    void initStructures();
    ~BamIndexer();

    template<typename Sink> std::streamsize write(Sink &snk, const char* s, std::streamsize n);
    template<typename Sink> bool flush(Sink& snk);
    void close();

private:
    void parseBgzfStream( const char* s, const std::streamsize src_size );
    void ProcessBgzfBlock();

    void ParseDecompressedBam();
    void addToBinIndex( const uint32_t bin, const VirtualOffset& virtualOffset, const VirtualOffset& virtualEndOffset );
    void addToLinearIndex( const uint32_t pos, const VirtualOffset& virtualOffset );

    void outputBaiHeader();
    void outputBaiFooter();
    void outputBaiChromosomeIndex();

private:
    Device baiSink_;
    bios::gzip_decompressor decompressor_;

    // BGZF parser
    enum {BGZF_STAGE_INIT, BGZF_STAGE_HEADER, BGZF_STAGE_BODY, BGZF_STAGE_FOOTER} bgzfParserStage_;
    uint32_t bgzfParserBytesNeeded_;
    std::vector<unsigned char> bgzfBuf_;
    uint64_t bgzfBlockCompressedOffset_;
    uint32_t uncompressedOffsetInBgzfBlock_;

    // From BGZF to BAM parsers
    std::vector<char> decompressedBam_;

    // BAM parser
    enum {
        BAM_STAGE_INIT
        , BAM_STAGE_HEADER
        , BAM_STAGE_SAM_HEADER_TEXT
        , BAM_STAGE_REF_SEQ_NUM
        , BAM_STAGE_REF_NAME_LENGTH
        , BAM_STAGE_REF_SEQ_INFO
        , BAM_STAGE_ALIGNMENT_BLOCK_SIZE
        , BAM_STAGE_ALIGNMENT_DATA
    } bamParserStage_;
    uint32_t bamParserBytesNeeded_;
    uint32_t bamParserStageLoopLeft_;
    VirtualOffset bamParserCurrentVirtualOffset_;
    VirtualOffset bamParserCurrentVirtualEndOffset_;
    VirtualOffset bamParserNextVirtualOffset_;
    uint64_t bamStatsMapped_, bamStatsNmapped_;
    int bamRefCount_;
    int lastProcessedRefId_;

    // Bin index
    uint32_t lastIndexedBin_;
    std::vector< std::vector< VirtualOffsetPair > > binIndex_;

    // Linear index
    std::vector< VirtualOffset > linearIndex_;
};


template<typename Device>
BamIndexer<Device>::BamIndexer(Device& baiSink)
    : baiSink_(baiSink)
    , decompressor_( bios::zlib::default_window_bits, GZIP_INTERNAL_RAM )
    , bgzfParserStage_( BGZF_STAGE_INIT )
    , bgzfParserBytesNeeded_( 0 )
    , bgzfBlockCompressedOffset_( 0 )
    , uncompressedOffsetInBgzfBlock_( 0 )
    , bamParserStage_( BAM_STAGE_INIT )
    , bamParserBytesNeeded_( 0 )
    , bamStatsMapped_( 0 )
    , bamStatsNmapped_( 0 )
    , lastProcessedRefId_( -1 )
    , lastIndexedBin_( 0 )
{
    binIndex_.resize( BAM_MAX_BIN );
    initStructures();
}

template<typename Device>
BamIndexer<Device>::BamIndexer(const BamIndexer& that)
    : baiSink_( that.baiSink_ )
    , decompressor_( that.decompressor_ )
    , bgzfParserStage_( that.bgzfParserStage_ )
    , bgzfParserBytesNeeded_( that.bgzfParserBytesNeeded_ )
    , bgzfBuf_( that.bgzfBuf_ )
    , bgzfBlockCompressedOffset_( that.bgzfBlockCompressedOffset_ )
    , uncompressedOffsetInBgzfBlock_( that.uncompressedOffsetInBgzfBlock_ )
    , bamParserStage_( that.bamParserStage_ )
    , bamParserBytesNeeded_( that.bamParserBytesNeeded_ )
    , bamParserCurrentVirtualOffset_( that.bamParserCurrentVirtualOffset_ )
    , bamParserCurrentVirtualEndOffset_( that.bamParserCurrentVirtualEndOffset_ )
    , bamParserNextVirtualOffset_( that.bamParserNextVirtualOffset_ )
    , bamStatsMapped_( that.bamStatsMapped_ )
    , bamStatsNmapped_( that.bamStatsNmapped_ )
    , lastProcessedRefId_( that.lastProcessedRefId_ )
    , lastIndexedBin_( that.lastIndexedBin_ )
    , binIndex_( that.binIndex_  )
    , linearIndex_( that.linearIndex_  )
{
    initStructures();
}


template<typename Device>
void BamIndexer<Device>::initStructures()
{
    ISAAC_ASSERT_MSG( binIndex_.size() == BAM_MAX_BIN, "Unexpected number of bins in Bam index" );
    for (uint32_t i=0; i<BAM_MAX_BIN; ++i)
    {
        binIndex_[i].reserve( MAX_CLUSTER_PER_INDEX_BIN );
    }

    linearIndex_.reserve( BAM_MAX_CONTIG_LENGTH / 16384 );
    decompressedBam_.reserve( MAX_UNCOMPRESSED_SIZE );
    bgzfBuf_.reserve( MAX_COMPRESSED_SIZE );
}


template<typename Device>
BamIndexer<Device>::~BamIndexer()
{
    if (bgzfBlockCompressedOffset_)
    {
        while (lastProcessedRefId_ != bamRefCount_)
        {
            ISAAC_ASSERT_MSG (lastProcessedRefId_ < bamRefCount_,
                              "Bam indexer processed more chromosomes than was declared in Bam header" );
            outputBaiChromosomeIndex();
            lastProcessedRefId_++;
        }
        outputBaiFooter();
    }
}


template<typename Device>
template<typename Sink>
std::streamsize BamIndexer<Device>::write(Sink &snk, const char* s, std::streamsize src_size)
{
    std::streamsize bytesWritten = bios::write(snk, s, src_size);
    ISAAC_ASSERT_MSG( bytesWritten == src_size, "Could not transfer all bytes from BAM source to BAM output" );
    parseBgzfStream( s, src_size );
    return bytesWritten;
}


template<typename Device>
void BamIndexer<Device>::parseBgzfStream( const char* inputBlock, const std::streamsize src_size )
{
    // Parse incoming compressed stream
    std::streamsize bytesProcessed( 0 );
    std::streamsize bytesLeft     ( src_size );

    while (bytesLeft >= bgzfParserBytesNeeded_)
    {
        ISAAC_ASSERT_MSG( (bgzfParserBytesNeeded_ > 0) || (bgzfParserStage_ == BGZF_STAGE_INIT),
                          "BGZF parser shouldn't be waiting for 0 bytes" );
        bgzfBuf_.insert( bgzfBuf_.end(), inputBlock + bytesProcessed, inputBlock + bytesProcessed + bgzfParserBytesNeeded_ );
        bytesProcessed += bgzfParserBytesNeeded_;
        bytesLeft -= bgzfParserBytesNeeded_;
        bgzfParserBytesNeeded_ = 0;

        switch (bgzfParserStage_)
        {
        case BGZF_STAGE_INIT:
        {
            ISAAC_ASSERT_MSG( bgzfBuf_.size() == 0, "BGZF parser was not initialised correctly" );
            bgzfParserBytesNeeded_ = 18;
            bgzfParserStage_ = BGZF_STAGE_HEADER;
            break;
        }
        case BGZF_STAGE_HEADER:
        {
            ISAAC_ASSERT_MSG( bgzfBuf_.size() == 18, "BGZF header was expected to be 18 bytes long but is different" );

            const uint32_t xLen  = (bgzfBuf_[11] << 8) + bgzfBuf_[10];
            const uint32_t bSize = (bgzfBuf_[17] << 8) + bgzfBuf_[16];

            bgzfParserBytesNeeded_ = bSize - xLen - 19;
            bgzfParserStage_ = BGZF_STAGE_BODY;
            break;
        }
        case BGZF_STAGE_BODY:
        {
            bgzfParserBytesNeeded_ = 8;
            bgzfParserStage_ = BGZF_STAGE_FOOTER;
            break;
        }
        case BGZF_STAGE_FOOTER:
        {
            uint32_t uncompressedSize = *(reinterpret_cast<uint32_t*>(&bgzfBuf_[4]));
            uncompressedSize++; // Avoid 'unused variable' warning in case warnings are not activated

            ProcessBgzfBlock();

            bgzfParserBytesNeeded_ = 18;
            bgzfBuf_.clear();
            bgzfParserStage_ = BGZF_STAGE_HEADER;
            break;
        }
        }
    }

    bgzfBuf_.insert( bgzfBuf_.end(), inputBlock + bytesProcessed, inputBlock + bytesProcessed + bytesLeft );
    bgzfParserBytesNeeded_ -= bytesLeft;
}


template<typename Device>
void BamIndexer<Device>::ProcessBgzfBlock()
{
    uint32_t lastDecompressedBamSize = decompressedBam_.size();
    bios::back_insert_device< std::vector<char> > decompressorSnk(decompressedBam_);
    bgzfBuf_[3] = '\0'; // reset FLG field of BGZF to discard extra subfields, then skip the extra subfield at bytes 11-17
    decompressor_.write(decompressorSnk, reinterpret_cast<const char*>(&bgzfBuf_[0]) , 10);
    decompressor_.write(decompressorSnk, reinterpret_cast<const char*>(&bgzfBuf_[18]), bgzfBuf_.size()-18);

    uint32_t bgzfCompressedSize = bgzfBuf_.size();
    uint32_t bgzfDecompressedSize = decompressedBam_.size() - lastDecompressedBamSize;
    bgzfDecompressedSize++; // Avoid 'unused variable' warning in case warnings are not activated

    ParseDecompressedBam();

    bgzfBlockCompressedOffset_ += bgzfCompressedSize;
}


template<typename Device>
void BamIndexer<Device>::ParseDecompressedBam()
{
    char *bamPtr = &decompressedBam_[0];
    std::streamsize bytesLeft( decompressedBam_.size() );

    while (bytesLeft >= bamParserBytesNeeded_)
    {
        uint32_t bytesToParse = bamParserBytesNeeded_;
        bamParserBytesNeeded_ = 0;

        switch (bamParserStage_)
        {
        case BAM_STAGE_INIT:
        {
            bamParserBytesNeeded_ = 8;
            bamParserStage_ = BAM_STAGE_HEADER;
            break;
        }
        case BAM_STAGE_HEADER:
        {
            ISAAC_ASSERT_MSG( bamPtr[0] == 'B', "Corrupted uncompressed BAM header" );
            ISAAC_ASSERT_MSG( bamPtr[1] == 'A', "Corrupted uncompressed BAM header" );
            ISAAC_ASSERT_MSG( bamPtr[2] == 'M', "Corrupted uncompressed BAM header" );
            ISAAC_ASSERT_MSG( bamPtr[3] == '\1', "Corrupted uncompressed BAM header" );

            const uint32_t lText = *reinterpret_cast< uint32_t* >(&bamPtr[4]);
            bamParserBytesNeeded_ = lText;
            bamParserStage_ = BAM_STAGE_SAM_HEADER_TEXT;
            break;
        }
        case BAM_STAGE_SAM_HEADER_TEXT:
        {
            bamParserBytesNeeded_ = 4;
            bamParserStage_ = BAM_STAGE_REF_SEQ_NUM;
            break;
        }
        case BAM_STAGE_REF_SEQ_NUM:
        {
            bamRefCount_ = *reinterpret_cast< uint32_t* >(&bamPtr[0]);

            ISAAC_ASSERT_MSG (bamRefCount_ > 0, "Invalid number of sequences in uncompressed BAM (n_ref=0)" );
            bamParserStageLoopLeft_ = bamRefCount_;
            bamParserBytesNeeded_ = 4;
            bamParserStage_ = BAM_STAGE_REF_NAME_LENGTH;
            break;
        }
        case BAM_STAGE_REF_NAME_LENGTH:
        {
            const uint32_t lName = *reinterpret_cast< uint32_t* >(&bamPtr[0]);

            ISAAC_ASSERT_MSG (lName > 0, "Invalid chromosome name in uncompressed BAM" );
            bamParserBytesNeeded_ = lName + 4;
            bamParserStage_ = BAM_STAGE_REF_SEQ_INFO;
            break;
        }
        case BAM_STAGE_REF_SEQ_INFO:
        {
            if (--bamParserStageLoopLeft_ > 0)
            {
                bamParserBytesNeeded_ = 4;
                bamParserStage_ = BAM_STAGE_REF_NAME_LENGTH;
            }
            else
            {
                outputBaiHeader();

                bamParserBytesNeeded_ = 4;
                bamParserStage_ = BAM_STAGE_ALIGNMENT_BLOCK_SIZE;
            }
            break;
        }
        case BAM_STAGE_ALIGNMENT_BLOCK_SIZE:
        {
            const uint32_t blockSize = *reinterpret_cast< uint32_t* >(&bamPtr[0]);

            if (bamParserNextVirtualOffset_.get())
            {
                bamParserCurrentVirtualOffset_ = bamParserNextVirtualOffset_;
            }
            else
            {
                bamParserCurrentVirtualOffset_.set( bgzfBlockCompressedOffset_, uncompressedOffsetInBgzfBlock_);
            }
            bamParserBytesNeeded_ = blockSize;
            bamParserStage_ = BAM_STAGE_ALIGNMENT_DATA;
            break;
        }
        case BAM_STAGE_ALIGNMENT_DATA:
        {
            struct BamAlignment
            {
                int refId;
                uint32_t pos;
                uint32_t binMqNl;
                uint32_t flagNc;
                uint32_t lSeq;
            } __attribute__ ((packed));

            const BamAlignment &alignment(*reinterpret_cast<BamAlignment*>(&bamPtr[0]));
            const uint32_t bin = alignment.binMqNl >> 16;
            const uint32_t flag = alignment.flagNc >> 16;

            while (alignment.refId != lastProcessedRefId_)
            {
                ISAAC_ASSERT_MSG (alignment.refId > lastProcessedRefId_,
                                  "Chromosome number in BAM alignment is greater than the number of chromosomes declared in BAM header" );
                if (-1 != lastProcessedRefId_)
                {
                    outputBaiChromosomeIndex();
                }
                lastProcessedRefId_++;
            }

            bamParserCurrentVirtualEndOffset_.set( bgzfBlockCompressedOffset_, uncompressedOffsetInBgzfBlock_ + bytesToParse);
            if ( bytesLeft == bytesToParse )
            {
                bamParserCurrentVirtualEndOffset_.set( bgzfBlockCompressedOffset_ + bgzfBuf_.size(), 0);
            }

            if (flag & BAM_FUNMAP)
            {
                // Update bamStats for samtools' special bin
                ++bamStatsNmapped_;
            } else {
                addToBinIndex( bin, bamParserCurrentVirtualOffset_, bamParserCurrentVirtualEndOffset_ );
                addToLinearIndex( alignment.pos, bamParserCurrentVirtualOffset_ );
                addToLinearIndex( alignment.pos + alignment.lSeq - 1, bamParserCurrentVirtualOffset_ );
                // Update bamStats for samtools' special bin
                ++bamStatsMapped_;
            }

            bamParserBytesNeeded_ = 4;
            bamParserStage_ = BAM_STAGE_ALIGNMENT_BLOCK_SIZE;
            break;
        }
        }

        uncompressedOffsetInBgzfBlock_ += bytesToParse;
        bamPtr    += bytesToParse;
        bytesLeft -= bytesToParse;
        if (bytesLeft > 0)
        {
            bamParserNextVirtualOffset_.set( bgzfBlockCompressedOffset_, uncompressedOffsetInBgzfBlock_);
        }
    }
    if (bytesLeft == 0)
    {
        bamParserNextVirtualOffset_.set( 0 );
    }

    decompressedBam_.erase( decompressedBam_.begin(), decompressedBam_.begin()+(bamPtr-&decompressedBam_[0]) );
    uncompressedOffsetInBgzfBlock_ = -decompressedBam_.size();
}


template<typename Device>
void BamIndexer<Device>::addToBinIndex( const uint32_t bin, const VirtualOffset& virtualOffset, const VirtualOffset& virtualEndOffset )
{
    ISAAC_ASSERT_MSG( binIndex_.size() == BAM_MAX_BIN, "Unexpected number of bins in Bam index" );
    ISAAC_ASSERT_MSG( bin < BAM_MAX_BIN, "Invalid bin number in uncompressed BAM" );

    if (binIndex_[bin].empty() || 
        ( bin != lastIndexedBin_ && binIndex_[bin].rbegin()->second.compressedOffset() != virtualOffset.compressedOffset() ))
    {
        binIndex_[bin].push_back( std::make_pair( virtualOffset, virtualEndOffset ) );
        lastIndexedBin_ = bin;
    }
    else
    {
        binIndex_[bin].rbegin()->second = virtualEndOffset;
    }
}


template<typename Device>
void BamIndexer<Device>::addToLinearIndex( const uint32_t pos, const VirtualOffset& virtualOffset )
{
    ISAAC_ASSERT_MSG( pos < BAM_MAX_CONTIG_LENGTH, "Alignment position greater than the maximum allowed by BAM index" );
    const uint32_t linearBin = pos>>14;
    if ( linearIndex_.size() <= linearBin )
    {
        const VirtualOffset& lastValue = linearIndex_.empty()?VirtualOffset():linearIndex_.back();
        while ( linearIndex_.size() <= linearBin )
        {
            linearIndex_.push_back( lastValue );
        }
        linearIndex_[linearBin] = virtualOffset;
    }
}


template<typename Device>
void BamIndexer<Device>::outputBaiHeader()
{
    bios::write(baiSink_, "BAI\1", 4);
    bios::write(baiSink_, reinterpret_cast<const char*>(&bamRefCount_), 4);
}


template<typename Device>
void BamIndexer<Device>::outputBaiFooter()
{
    // output number of coor-less reads (special samtools field)
    const uint64_t zero = 0;
    bios::write(baiSink_, reinterpret_cast<const char*>(&zero), 8);
}


template<typename Device>
void BamIndexer<Device>::outputBaiChromosomeIndex()
{
    struct {
        uint32_t binNum;
        uint32_t nClusters;
        uint64_t offBeg, offEnd;
        uint64_t mapped, nmapped;
    } __attribute__ ((packed))
          specialBin = { BAM_MAX_BIN, 2, -1U, 0, bamStatsMapped_, bamStatsNmapped_ };

    const uint32_t nBin = std::count_if( binIndex_.begin(), 
                                             binIndex_.end(),
                                             boost::bind(&std::vector< VirtualOffsetPair >::empty, _1) == false )
        + 1; // Add samtools' special bin to the count
    bios::write(baiSink_, reinterpret_cast<const char*>(&nBin), 4);

    uint32_t i=0;
    BOOST_FOREACH( const std::vector< VirtualOffsetPair >& binIndexEntry, binIndex_ )
    {
        if ( !binIndexEntry.empty() )
        {
            bios::write(baiSink_, reinterpret_cast<const char*>(&i), 4);
            const uint32_t nChunk = binIndexEntry.size();
            bios::write(baiSink_, reinterpret_cast<const char*>(&nChunk), 4);
            bios::write(baiSink_, reinterpret_cast<const char*>(&binIndexEntry[0]), nChunk*16);

            // Fill in samtools' "specialBin" bamStats
            if (specialBin.offBeg > binIndexEntry[0].first.get())
            {
                specialBin.offBeg = binIndexEntry[0].first.get();
            }
            if (specialBin.offEnd < binIndexEntry[nChunk-1].second.get())
            {
                specialBin.offEnd = binIndexEntry[nChunk-1].second.get();
            }
        }
        ++i;
    }

    // Writing special samtools bin
    bios::write(baiSink_, reinterpret_cast<const char*>(&specialBin), sizeof(specialBin));

    // n_intv
    const uint32_t nIntv = linearIndex_.size();
    bios::write(baiSink_, reinterpret_cast<const char*>(&nIntv), 4);
    bios::write(baiSink_, reinterpret_cast<const char*>(&linearIndex_[0]), nIntv*8);


    // reset variables to make them ready to process the next chromosome
    bamStatsMapped_ = bamStatsNmapped_ = 0;
    ISAAC_ASSERT_MSG( binIndex_.size() == BAM_MAX_BIN, "Unexpected number of bins in Bam index" );
    BOOST_FOREACH( std::vector< VirtualOffsetPair >& binIndexEntry, binIndex_)
    {
        binIndexEntry.clear();
    }
    linearIndex_.clear();
}


template<typename Device>
void BamIndexer<Device>::close()
{
    bios::close( baiSink_ );
}


template<typename Device>
template<typename Sink>
bool BamIndexer<Device>::flush(Sink& snk)
{
    bool r1 = bios::flush( snk );
    bool r2 = bios::flush( baiSink_ );
    return r1 && r2;
}





////////////////////////////////

typedef std::pair<uint32_t,uint32_t> Chunk;

struct UnresolvedBinIndexChunk
{
    UnresolvedBinIndexChunk() {}
    UnresolvedBinIndexChunk(UnresolvedOffset startPos,
                            UnresolvedOffset endPos,
                            uint32_t  bin,
                            uint32_t  refId)
        : startPos(startPos)
        , endPos  (endPos)
        , bin     (bin)
        , refId   (refId)
        {}

    UnresolvedOffset startPos;
    UnresolvedOffset endPos;
    uint32_t  bin;
    uint32_t  refId;
};

class BamIndexPart
{
    static const uint32_t BAM_INDEXER_MAX_CHUNKS = BAM_MAX_BIN * MAX_CLUSTER_PER_INDEX_BIN;         
    static const uint32_t BAM_MIN_CHUNK_GAP = 32768;

public:
    BamIndexPart();
    void processFragment( const build::FragmentAccessorBamAdapter& alignment, uint32_t serializedLength );

//private:
    void initStructures();
    void addToBinIndexChunks( const UnresolvedOffset virtualOffset, const UnresolvedOffset virtualEndOffset, const uint32_t bin, const uint32_t refId );
    void addToLinearIndex( const uint32_t pos, const UnresolvedOffset virtualOffset );


    UnresolvedOffset localUncompressedOffset_;

    // Bin index
    std::vector<UnresolvedBinIndexChunk> chunks_;

    // Linear index
    std::vector< UnresolvedOffset > linearIndex_;

    // Stats reported in last bin
    uint64_t bamStatsMapped_, bamStatsNmapped_;
};


class BamIndex
{
public:
    // Creates invalid object which is not to be used
    BamIndex();
    // Creates proper object with output file attached
    BamIndex(const boost::filesystem::path &bamPath, const uint32_t bamRefCount, const uint32_t bamHeaderCompressedLength);
    void processIndexPart(const bam::BamIndexPart &bamIndexPart,
                          const std::vector<char> &bgzfBuffer);

    void flush()
    {
        outputIndexFile();
    }

private:
    void initStructures();
    void outputIndexFile();
    void outputBaiHeader();
    void outputBaiFooter();
    void outputBaiChromosomeIndex();

    void printBgzfInfo( const std::vector<char> bgzfBuffer );
    void printBamIndexPartInfo( const bam::BamIndexPart &bamIndexPart );

    void mergeBinIndex( const std::vector<UnresolvedBinIndexChunk>& binIndexChunks, const std::vector<char> &bgzfBuffer );
    void mergeLinearIndex( const std::vector<UnresolvedOffset>& linearIndexToMerge, const std::vector<char> &bgzfBuffer );
    void addToBinIndex( const UnresolvedBinIndexChunk& chunk, const std::vector<char> &bgzfBuffer );
    void clearStructures();
    void resetBgzfParsing();
    VirtualOffset resolveOffset(UnresolvedOffset unresolvedPos, const std::vector<char> &bgzfBuffer);

    uint32_t bamRefCount_;
    uint32_t lastProcessedRefId_;
    std::ofstream baiStream_;

    // Bin index
    std::vector< std::vector< VirtualOffsetPair > > binIndex_;
    // As binIndex_ might have quite a few entries, speed up cleanup and counting by tracking whether
    // anything has been put into it with binIndexEmpty_
    bool binIndexEmpty_;

    // Linear index
    std::vector< VirtualOffset > linearIndex_;

    // Stats reported in last bin
    uint64_t bamStatsMapped_, bamStatsNmapped_, bamStatsGlobalNoCoordinates_;

    uint64_t positionInBam_;
    uint64_t currentBgzfBlockCompressedPosition_;
    uint64_t currentBgzfBlockUncompressedPosition_;
    uint32_t currentBgzfBlockCompressedSize_;
    uint32_t currentBgzfBlockUncompressedSize_;
};


} // namespace bam
} // namespace isaac


#endif // iSAAC_BAM_BAM_INDEXER_HH
