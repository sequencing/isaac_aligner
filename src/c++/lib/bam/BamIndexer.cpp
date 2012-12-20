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
 ** \file BamIndexer.cpp
 **
 ** \brief Implementation of BAM indexing
 **
 ** \author Lilian Janin
 **/

#include "bam/BamIndexer.hh"
#include "alignment/Cigar.hh"


namespace isaac
{
namespace bam
{


BamIndexPart::BamIndexPart()
    : localUncompressedOffset_( 0 )
    , bamStatsMapped_( 0 )
    , bamStatsNmapped_( 0 )
{
    initStructures();
}

void BamIndexPart::initStructures()
{
    chunks_.reserve( BAM_INDEXER_MAX_CHUNKS );
    linearIndex_.reserve( BAM_MAX_CONTIG_LENGTH / 16384 );
}

void BamIndexPart::processFragment( const build::FragmentAccessorBamAdapter& alignment, uint32_t serializedLength )
{
    const uint32_t flag(alignment.flag());

    if (alignment.pos() >= 0)
    {
        uint32_t observedLength = alignment.observedLength();
        const uint32_t bin(bam_reg2bin(alignment.pos(), alignment.pos() + alignment.seqLen())); // it would be more correct to use observedLength instead of alignment.seqLen(), but samtools is doing it this way.

        addToBinIndexChunks( localUncompressedOffset_, localUncompressedOffset_ + serializedLength, bin, alignment.refId() );
        addToLinearIndex( alignment.pos(), localUncompressedOffset_ );
        if (observedLength > 0)
        {
            addToLinearIndex( alignment.pos() + observedLength - 1, localUncompressedOffset_ );
        }
    }

    // Update bamStats for samtools' special bin
    if (flag & BAM_FUNMAP)
    {
        ++bamStatsNmapped_;
    }
    else
    {
        ++bamStatsMapped_;
    }

    localUncompressedOffset_ += serializedLength;
}

void BamIndexPart::addToBinIndexChunks( const UnresolvedOffset virtualOffset, const UnresolvedOffset virtualEndOffset, const uint32_t bin, const uint32_t refId )
{
    ISAAC_ASSERT_MSG( bin < BAM_MAX_BIN, "Invalid bin number in uncompressed BAM" );

    if (!chunks_.empty() &&
        bin == chunks_.back().bin &&
        refId == chunks_.back().refId)
    {
        chunks_.back().endPos = virtualEndOffset;
    }
    else if (chunks_.size() >= 2 &&
             bin == chunks_[chunks_.size()-2].bin &&
             refId == chunks_[chunks_.size()-2].refId &&
             (chunks_[chunks_.size()-2].endPos + BAM_MIN_CHUNK_GAP) > virtualEndOffset)
    {
        // Chunk reduction around the boundary of two adjacent bins
        chunks_[chunks_.size()-2].endPos = virtualEndOffset;
    }
    else
    {
        chunks_.push_back(UnresolvedBinIndexChunk(virtualOffset, virtualEndOffset, bin, refId));
    }
}

void BamIndexPart::addToLinearIndex( const uint32_t pos, const UnresolvedOffset virtualOffset )
{
    if (pos >= BAM_MAX_CONTIG_LENGTH)
    {
        ISAAC_ASSERT_MSG( pos < BAM_MAX_CONTIG_LENGTH, "Alignment position greater than the maximum allowed by BAM index: " << pos);
    }
    const uint32_t linearBin = pos>>14;
    if ( linearIndex_.size() <= linearBin )
    {
        const UnresolvedOffset lastValue = linearIndex_.empty()?0xFFFFFFFFFFFFFFFF:linearIndex_.back();
        while ( linearIndex_.size() <= linearBin )
        {
            linearIndex_.push_back( lastValue );
        }
        linearIndex_[linearBin] = virtualOffset;
    }
}


BamIndex::BamIndex()
    : bamRefCount_( 0 )
    , lastProcessedRefId_( 0xFFFFFFFF )
    , baiStream_()
    , binIndex_( BAM_MAX_BIN )
    , bamStatsMapped_( 0 )
    , bamStatsNmapped_( 0 )
    , positionInBam_( 0 )
    , currentBgzfBlockCompressedPosition_( 0 )
    , currentBgzfBlockUncompressedPosition_( 0 )
    , currentBgzfBlockCompressedSize_( 0 )
    , currentBgzfBlockUncompressedSize_( 0 )
{
}


BamIndex::BamIndex(const boost::filesystem::path &bamPath, const uint32_t bamRefCount, const uint32_t bamHeaderCompressedLength)
    : bamRefCount_( bamRefCount )
    , lastProcessedRefId_( 0xFFFFFFFF )
    , baiStream_( (bamPath.string() + ".bai").c_str() )
    , binIndex_( BAM_MAX_BIN )
    , bamStatsMapped_( 0 )
    , bamStatsNmapped_( 0 )
    , bamStatsGlobalNoCoordinates_( 0 )
    , positionInBam_( bamHeaderCompressedLength )
    , currentBgzfBlockCompressedPosition_( 0 )
    , currentBgzfBlockUncompressedPosition_( 0 )
    , currentBgzfBlockCompressedSize_( 0 )
    , currentBgzfBlockUncompressedSize_( 0 )
{
    ISAAC_ASSERT_MSG( baiStream_.good(), "Error opening bam index file for writing" );
    initStructures();
    outputBaiHeader();
}

BamIndex::~BamIndex()
{
    outputIndexFile();
}

void BamIndex::initStructures()
{
    ISAAC_ASSERT_MSG( binIndex_.size() == BAM_MAX_BIN, "Unexpected number of bins in Bam index" );
    for (uint32_t i=0; i<BAM_MAX_BIN; ++i)
    {
        binIndex_[i].reserve( MAX_CLUSTER_PER_INDEX_BIN );
    }
    linearIndex_.reserve( BAM_MAX_CONTIG_LENGTH / 16384 );
}

void BamIndex::outputIndexFile()
{
    if (lastProcessedRefId_ == 0xFFFFFFFF)
    {
        lastProcessedRefId_ = 0;
    }

    while (lastProcessedRefId_ != bamRefCount_)
    {
        ISAAC_ASSERT_MSG (lastProcessedRefId_ < bamRefCount_,
                          "Bam indexer processed more chromosomes than was declared in Bam header" );
        outputBaiChromosomeIndex();
        lastProcessedRefId_++;
    }
    outputBaiFooter();
}

void BamIndex::outputBaiHeader()
{
    baiStream_.write("BAI\1", 4);
    baiStream_.write(reinterpret_cast<const char*>(&bamRefCount_), 4);
}

void BamIndex::outputBaiChromosomeIndex()
{
    struct {
        uint32_t binNum;
        uint32_t nClusters;
        uint64_t offBeg, offEnd;
        uint64_t mapped, nmapped;
    } __attribute__ ((packed))
          specialBin = { BAM_MAX_BIN, 2, 0, 0, bamStatsMapped_, bamStatsNmapped_ };

    uint32_t nBin = std::count_if( binIndex_.begin(), 
                                   binIndex_.end(),
                                   boost::bind(&std::vector< VirtualOffsetPair >::empty, _1) == false );

    if (nBin > 0 || bamStatsMapped_ > 0 || bamStatsNmapped_ > 0)
    {
        ++nBin; // Add samtools' special bin to the count
        baiStream_.write(reinterpret_cast<const char*>(&nBin), 4);

        unsigned int i=0;
        BOOST_FOREACH( const std::vector< VirtualOffsetPair >& binIndexEntry, binIndex_ )
        {
            if ( !binIndexEntry.empty() )
            {
                baiStream_.write(reinterpret_cast<const char*>(&i), 4);
                const uint32_t nChunk = binIndexEntry.size();
                baiStream_.write(reinterpret_cast<const char*>(&nChunk), 4);
                baiStream_.write(reinterpret_cast<const char*>(&binIndexEntry[0]), nChunk*16);

                // Fill in samtools' "specialBin" bamStats
                if (specialBin.offBeg > binIndexEntry[0].first.get() || specialBin.offBeg == 0)
                {
                    specialBin.offBeg = binIndexEntry[0].first.get();
                }
                if (specialBin.offEnd < binIndexEntry[nChunk-1].second.get() || specialBin.offEnd == 0)
                {
                    specialBin.offEnd = binIndexEntry[nChunk-1].second.get();
                }
            }
            ++i;
        }

        // Write special samtools bin
        baiStream_.write(reinterpret_cast<const char*>(&specialBin), sizeof(specialBin));
    }
    else
    {
        baiStream_.write(reinterpret_cast<const char*>(&nBin), 4); // nBin==0
    }

    // Write linear index
    const uint32_t nIntv = linearIndex_.size();
    baiStream_.write(reinterpret_cast<const char*>(&nIntv), 4);
    baiStream_.write(reinterpret_cast<const char*>(&linearIndex_[0]), nIntv*8);

    // reset variables to make them ready to process the next chromosome
    clearStructures();
}

void BamIndex::outputBaiFooter()
{
    // output number of coor-less reads (special samtools field)
    baiStream_.write(reinterpret_cast<const char*>(&bamStatsGlobalNoCoordinates_), 8);
}

void BamIndex::processIndexPart(const bam::BamIndexPart &bamIndexPart,
                                const std::vector<char> &bgzfBuffer)
{
    if (bgzfBuffer.empty())
    {
        return;
    }

    if (!bamIndexPart.chunks_.empty())
    {
        // Block of mapped reads
        uint32_t refId = bamIndexPart.chunks_[0].refId;
        while (lastProcessedRefId_ != refId)
        {
            if (lastProcessedRefId_ == 0xFFFFFFFF)
            {
                clearStructures();
                lastProcessedRefId_ = 0;
            }
            else
            {
                ISAAC_ASSERT_MSG (lastProcessedRefId_ < refId,
                                  "Bam indexer tries to process more chromosomes than was declared in Bam header" );
                outputBaiChromosomeIndex();
                lastProcessedRefId_++;
            }
        }

        resetBgzfParsing();
        mergeBinIndex( bamIndexPart.chunks_, bgzfBuffer );
        mergeLinearIndex( bamIndexPart.linearIndex_, bgzfBuffer );

        bamStatsMapped_ += bamIndexPart.bamStatsMapped_;
        bamStatsNmapped_ += bamIndexPart.bamStatsNmapped_;
    }
    else
    {
        // Block of unmapped reads
        bamStatsGlobalNoCoordinates_ += bamIndexPart.bamStatsNmapped_;
    }

    // Add offset for next index part
    positionInBam_ += bgzfBuffer.size();
}

void BamIndex::printBgzfInfo( const std::vector<char> bgzfBuffer )
{
    uint64_t compressedSize = bgzfBuffer.size();
    uint64_t uncompressedSize = 0;
    baiStream_ << "total compressedSize:   " << compressedSize << std::endl;
    baiStream_ << "uncompressedSize: " << uncompressedSize << std::endl;

    uint64_t posInBgzf = 0;
    while (posInBgzf < bgzfBuffer.size())
    {
        baiStream_ << "** posInBgzf:   " << posInBgzf << std::endl;
        ISAAC_ASSERT_MSG( bgzfBuffer[posInBgzf+0] == '\x1f', "Error while parsing BGZF block during indexing: invalid byte 0" );
        ISAAC_ASSERT_MSG( bgzfBuffer[posInBgzf+1] == '\x8b', "Error while parsing BGZF block during indexing: invalid byte 1" );
        ISAAC_ASSERT_MSG( bgzfBuffer[posInBgzf+2] == '\x08', "Error while parsing BGZF block during indexing: invalid byte 2" );
        ISAAC_ASSERT_MSG( bgzfBuffer[posInBgzf+3] == '\x04', "Error while parsing BGZF block during indexing: invalid byte 3" );
        ISAAC_ASSERT_MSG( bgzfBuffer[posInBgzf+12] == '\x42', "Error while parsing BGZF block during indexing: invalid byte 12" );
        ISAAC_ASSERT_MSG( bgzfBuffer[posInBgzf+13] == '\x43', "Error while parsing BGZF block during indexing: invalid byte 13" );
        uint16_t compressedBlockSize   = *((uint16_t*)(&bgzfBuffer[posInBgzf+16]));
        uint32_t uncompressedBlockSize = *((uint32_t*)(&bgzfBuffer[posInBgzf+compressedBlockSize-3]));
        baiStream_ << "compressedBlockSize: "   << compressedBlockSize   << std::endl;
        baiStream_ << "uncompressedBlockSize: " << uncompressedBlockSize << std::endl;
        posInBgzf += compressedBlockSize + 1;
    }
    ISAAC_ASSERT_MSG( posInBgzf == bgzfBuffer.size(), "Error while parsing BGZF block during indexing: end of last block doesn't match end of buffer" );
}

void BamIndex::printBamIndexPartInfo( const bam::BamIndexPart &bamIndexPart )
{
    baiStream_ << "bamIndexPart.localUncompressedOffset_: " << bamIndexPart.localUncompressedOffset_ << std::endl;
    baiStream_ << "bamIndexPart.bamStatsMapped_: " << bamIndexPart.bamStatsMapped_ << std::endl;
    baiStream_ << "bamIndexPart.bamStatsNmapped_: " << bamIndexPart.bamStatsNmapped_ << std::endl;
    uint32_t binNum = 0;
    BOOST_FOREACH( const UnresolvedBinIndexChunk& chunk, bamIndexPart.chunks_ )
    {
        baiStream_ << "chunk[]: " << chunk.startPos << ", " << chunk.endPos << ", " << chunk.bin << ", " << chunk.refId << std::endl;
    }
    binNum = 0;
    BOOST_FOREACH( const UnresolvedOffset val, bamIndexPart.linearIndex_)
    {
        if (val)
        {
            baiStream_ << "linearIndex_[" << binNum << "]: " << val << std::endl;
        }
        ++binNum;
    }
}

void BamIndex::mergeBinIndex( const std::vector<UnresolvedBinIndexChunk>& binIndexChunks, const std::vector<char> &bgzfBuffer )
{
    BOOST_FOREACH( const UnresolvedBinIndexChunk& chunk, binIndexChunks )
    {
        addToBinIndex( chunk, bgzfBuffer );
    }
}

void BamIndex::mergeLinearIndex( const std::vector<UnresolvedOffset>& linearIndexToMerge, const std::vector<char> &bgzfBuffer )
{
    if (linearIndex_.size() < linearIndexToMerge.size())
    {
        linearIndex_.resize( linearIndexToMerge.size() );
    }
    for (unsigned i = 0; i < linearIndexToMerge.size(); ++i)
    {
        if (linearIndexToMerge[i] != 0xFFFFFFFFFFFFFFFF)
        {
            VirtualOffset off = resolveOffset( linearIndexToMerge[i], bgzfBuffer );
            if (off.get() < linearIndex_[i].get() || linearIndex_[i].get() == 0)
            {
                linearIndex_[i] = off;
            }
        }
    }
}

void BamIndex::addToBinIndex( const UnresolvedBinIndexChunk& chunk, const std::vector<char> &bgzfBuffer )
{
    ISAAC_ASSERT_MSG( binIndex_.size() == BAM_MAX_BIN, "Unexpected number of bins in Bam index" );
    ISAAC_ASSERT_MSG( chunk.bin < BAM_MAX_BIN, "Invalid bin number in uncompressed BAM" );

    VirtualOffset start = resolveOffset(chunk.startPos, bgzfBuffer);
    VirtualOffset end   = resolveOffset(chunk.endPos, bgzfBuffer);

    if (!binIndex_[chunk.bin].empty() && binIndex_[chunk.bin].back().second.compressedOffset() == start.compressedOffset())
    {
        // Small chunks reduction
        binIndex_[chunk.bin].back().second = end;
    }
    else
    {
        binIndex_[chunk.bin].push_back(std::make_pair(start, end));
    }
}

void BamIndex::clearStructures()
{
    bamStatsMapped_ = bamStatsNmapped_ = 0;
    ISAAC_ASSERT_MSG( binIndex_.size() == BAM_MAX_BIN, "Unexpected number of bins in Bam index" );
    BOOST_FOREACH( std::vector< VirtualOffsetPair >& binIndexEntry, binIndex_)
    {
        binIndexEntry.clear();
    }
    linearIndex_.clear();

    resetBgzfParsing();
}

void BamIndex::resetBgzfParsing()
{
    currentBgzfBlockCompressedPosition_ = 0;
    currentBgzfBlockUncompressedPosition_ = 0;
    currentBgzfBlockCompressedSize_ = 0;
    currentBgzfBlockUncompressedSize_ = 0;
}

VirtualOffset BamIndex::resolveOffset(UnresolvedOffset unresolvedPos, const std::vector<char> &bgzfBuffer)
{
    VirtualOffset result;
    if ( unresolvedPos < currentBgzfBlockUncompressedPosition_ )
    {
        // This could be made more efficient
        resetBgzfParsing();
    }
    while (unresolvedPos >= currentBgzfBlockUncompressedPosition_ + currentBgzfBlockUncompressedSize_)
    {
        currentBgzfBlockCompressedPosition_ += currentBgzfBlockCompressedSize_;
        currentBgzfBlockUncompressedPosition_ += currentBgzfBlockUncompressedSize_;
        if (currentBgzfBlockCompressedPosition_ == bgzfBuffer.size())
        {
            currentBgzfBlockCompressedSize_ = 0;
            currentBgzfBlockUncompressedSize_ = 0;
            break;
        }
        ISAAC_ASSERT_MSG( currentBgzfBlockCompressedPosition_+13 < bgzfBuffer.size(), "Error while parsing BGZF block during indexing: trying to read past end of buffer" );
        ISAAC_ASSERT_MSG( bgzfBuffer[currentBgzfBlockCompressedPosition_+0] == '\x1f', "Error while parsing BGZF block during indexing: invalid byte 0" );
        ISAAC_ASSERT_MSG( bgzfBuffer[currentBgzfBlockCompressedPosition_+1] == '\x8b', "Error while parsing BGZF block during indexing: invalid byte 1" );
        ISAAC_ASSERT_MSG( bgzfBuffer[currentBgzfBlockCompressedPosition_+2] == '\x08', "Error while parsing BGZF block during indexing: invalid byte 2" );
        ISAAC_ASSERT_MSG( bgzfBuffer[currentBgzfBlockCompressedPosition_+3] == '\x04', "Error while parsing BGZF block during indexing: invalid byte 3" );
        ISAAC_ASSERT_MSG( bgzfBuffer[currentBgzfBlockCompressedPosition_+12] == '\x42', "Error while parsing BGZF block during indexing: invalid byte 12" );
        ISAAC_ASSERT_MSG( bgzfBuffer[currentBgzfBlockCompressedPosition_+13] == '\x43', "Error while parsing BGZF block during indexing: invalid byte 13" );
        uint16_t compressedBlockSize   = *((uint16_t*)(&bgzfBuffer[currentBgzfBlockCompressedPosition_+16]));
        uint32_t uncompressedBlockSize = *((uint32_t*)(&bgzfBuffer[currentBgzfBlockCompressedPosition_+compressedBlockSize-3]));
        currentBgzfBlockCompressedSize_ = compressedBlockSize + 1;
        currentBgzfBlockUncompressedSize_ = uncompressedBlockSize;
    }

    result.set(currentBgzfBlockCompressedPosition_ + positionInBam_, unresolvedPos - currentBgzfBlockUncompressedPosition_);
    return result;
}


} // namespace bam
} // namespace isaac
