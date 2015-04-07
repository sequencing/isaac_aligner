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
 ** \file BamLoader.hh
 **
 ** Component to read Bam files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BAM_BAM_PARSER_HH
#define iSAAC_BAM_BAM_PARSER_HH

#include "common/Endianness.hh"
#include "common/Exceptions.hh"
#include "flowcell/ReadMetadata.hh"
#include "oligo/Nucleotides.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace bam
{

struct BamParserException : common::IoException
{
    BamParserException(const std::string &message) : common::IoException(EINVAL, message){}
};


inline char bamBase(const unsigned char bamSeq)
{
    static const unsigned char BAM_BASES[] = {'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
    const unsigned char bamBase = BAM_BASES[bamSeq];
    return bamBase;
}

inline unsigned char bamToBcl(const unsigned char qual, const unsigned char bamSeq)
{
    static const oligo::Translator translator = oligo::getTranslator();
    const unsigned char q = 0xFF == qual ? 0 : std::min<unsigned char>(qual, 0x3f);
    const unsigned char base = bamBase(bamSeq);
//            std::cerr << bamBase;
//            std::cerr << char(q + 33);
    const unsigned char baseValue = translator[base];
    return oligo::invalidOligo == baseValue ? 0 : (baseValue | (q << 2));
}

struct BamBlockHeader : boost::noncopyable
{
protected:
    int refID;
    int pos;
    unsigned bin_mq_nl;
    unsigned flag_nc;
    int l_seq;
    int next_refID;
    int next_pos;
    int tlen;
    char read_name[0];
public:

    static const unsigned MULTI_SEGMENT = 0x01 << 16;
    static const unsigned REV_COMPL = 0x10 << 16;
    static const unsigned FIRST_SEGMENT = 0x40 << 16;
    static const unsigned LAST_SEGMENT = 0x80 << 16;
    static const unsigned VERNDOR_FAILED = 0x200 << 16;
    static const unsigned SECONDARY_ALIGNMENT = 0x100 << 16;
    static const unsigned SUPPLEMENTARY_ALIGNMENT = 0x800 << 16;

    unsigned char getReadNameLength() const {return bin_mq_nl;}
    unsigned short getCigarLength() const {return flag_nc;}
    const char *nameBegin() const {return read_name ;}
    const char *nameEnd() const {return read_name + getReadNameLength();}
    const unsigned *getCigar() const {return reinterpret_cast<const unsigned*>(read_name + getReadNameLength());}
    const unsigned char *getSeq() const {return reinterpret_cast<const unsigned char *>(getCigar() + getCigarLength());}
    const unsigned char *getQual() const {return getSeq() + (getLSeq() + 1) / 2;}
    bool isPaired() const {return flag_nc & MULTI_SEGMENT;}
    bool isReverse() const {return flag_nc & REV_COMPL;}
    bool isReadOne() const
    {
        // FIRST_SEGMENT is not set for single-ended data
        return !(flag_nc & LAST_SEGMENT);
    }

    bool isPf() const {return !(flag_nc & VERNDOR_FAILED);}
    bool isSupplementaryAlignment() const {return flag_nc & SUPPLEMENTARY_ALIGNMENT;}
    bool isSecondaryAlignment() const {return flag_nc & SECONDARY_ALIGNMENT;}

    int getRefId() const {return common::extractLittleEndian<int>(&refID);}
    int getNextRefId() const {return common::extractLittleEndian<int>(&next_refID);}

    int getPos() const {return common::extractLittleEndian<int>(&pos);}
    int getNextPos() const {return common::extractLittleEndian<int>(&next_pos);}
    int getLSeq() const {return common::extractLittleEndian<int>(&l_seq);}

    friend std::ostream &seqToStream(std::ostream &os, const BamBlockHeader &bamBlockHeader)
    {
        const unsigned char * const pSeq = bamBlockHeader.getSeq();
        unsigned seqLen = bamBlockHeader.getLSeq();

        unsigned seqOffset = 0;
        while(seqLen--)
        {
            os << bamBase((*(pSeq + seqOffset / 2) >> (4 * ((seqOffset + 1) % 2))) & 0x0F);
            ++seqOffset;
        }

        return os;
    }

    friend std::ostream &operator << (std::ostream &os, const BamBlockHeader &bamBlockHeader)
    {
        return os << "BamBlockHeader("  << bamBlockHeader.read_name <<  "," <<
            bamBlockHeader.refID << ":" << bamBlockHeader.pos <<
            ")";
//        return seqToStream(os, bamBlockHeader) << ")";
    }
};

class BamParser
{
    unsigned headerBytesToSkip_;
    int referenceSequencesToSkip_;
public:
    BamParser() :
        headerBytesToSkip_(-1U),
        referenceSequencesToSkip_(-1) {}

    void reset()
    {
        headerBytesToSkip_ = -1U;
        referenceSequencesToSkip_ = -1;
    }

    template <typename CollectorT>
    bool parse(
        std::vector<char>::const_iterator &uncompressedIt,
        std::vector<char>::const_iterator uncompressedEnd,
        CollectorT &collector)
    {
        bool moreDataNeeded = true;
        if (headerBytesToSkip_)
        {
            skipHeader(uncompressedIt, uncompressedEnd);
        }

        if (!headerBytesToSkip_)
        {
            if (referenceSequencesToSkip_)
            {
                skipReferences(uncompressedIt, uncompressedEnd);
            }

            if (!referenceSequencesToSkip_)
            {
                while(uncompressedEnd != uncompressedIt)
                {
                    std::vector<char>::const_iterator last = uncompressedIt;
                    moreDataNeeded = parseBamRecord(uncompressedIt, uncompressedEnd, collector);
                    if(!moreDataNeeded || last == uncompressedIt)
                    {
                        // if the uncompressedIt did not move, means we can't parse from this point.
                        // if moreDataNeeded is false, the collector cannot accept anymore data
                        break;
                    }
                }
            }
        }

        return moreDataNeeded;
    }


private:
    bool skipHeader(
        std::vector<char>::const_iterator &it,
        std::vector<char>::const_iterator end);

    bool skipReferences(
        std::vector<char>::const_iterator &it,
        std::vector<char>::const_iterator end);

    template <typename ProcessorT>
    bool parseBamRecord(
        std::vector<char>::const_iterator &it,
        std::vector<char>::const_iterator end,
        ProcessorT &process)
    {
        static const unsigned BLOCK_SIZE_WIDTH = 4;
        const std::vector<char>::const_iterator blockIt = it;
        int block_size = 0;
        if (std::size_t(std::distance(it, end)) < BLOCK_SIZE_WIDTH)
        {
            return true;
        }

        it = common::extractLittleEndian(it, block_size);

        if (std::distance(it, end) < block_size)
        {
            it = blockIt;
            return true;
        }

        ISAAC_ASSERT_MSG(block_size >= int(sizeof(BamBlockHeader)), "bam record size is smaller than the minimum required block_size:" << block_size << " sizeof(BamBlockHeader):" << sizeof(BamBlockHeader));
        const BamBlockHeader &block = *reinterpret_cast<const BamBlockHeader *>(&*it);

        const bool lastBlock = std::size_t(std::distance(it + block_size, end)) <= BLOCK_SIZE_WIDTH ||
            std::distance(it + block_size + BLOCK_SIZE_WIDTH, end) < (common::extractLittleEndian<unsigned>(it + block_size));
        const bool ret = process(block, lastBlock);
        it += block_size;

        return ret;
    }

};

/**
 * \brief extracts up to corresponding read length bcl sequence. if bam sequence is shorter than read length, the
 *        rest is padded with 0
 */
template <typename InsertIt, typename UnaryFun>
InsertIt extractBcl(
    const BamBlockHeader &bamBlock,
    InsertIt insertIt,
    UnaryFun translate,
    const flowcell::ReadMetadata &readMetadata)
{
    const unsigned char * const pSeq = bamBlock.getSeq();
    const unsigned char *pQual = bamBlock.getQual();
    const unsigned char * const pQualEnd = pQual + bamBlock.getLSeq();
    unsigned currentCycle = readMetadata.getFirstReadCycle();
    std::vector<unsigned>::const_iterator cycleIterator = readMetadata.getCycles().begin();

    unsigned seqOffset = 0;
    while(pQualEnd != pQual && readMetadata.getCycles().end() != cycleIterator)
    {
        if (*cycleIterator == currentCycle)
        {
            *insertIt++ = translate(bamToBcl(*pQual,
                                             (*(pSeq + seqOffset / 2) >> (4 * ((seqOffset + 1) % 2))) & 0x0F
                                             ));
            ++cycleIterator;
        }
        ++pQual;
        ++seqOffset;
        ++currentCycle;
    }
    insertIt = std::fill_n(insertIt, std::distance(cycleIterator, readMetadata.getCycles().end()), 0);
    return insertIt;
}

template <typename InsertIt>
InsertIt extractForwardBcl(
    const BamBlockHeader &bamBlock,
    InsertIt insertIt,
    const flowcell::ReadMetadata &readMetadata)
{
    return extractBcl(bamBlock, insertIt, &boost::cref<unsigned char>, readMetadata);
}

template <typename RandomAccessIt>
RandomAccessIt extractReverseBcl(
    const BamBlockHeader &bamBlock,
    RandomAccessIt randomAccessIt,
    const flowcell::ReadMetadata &readMetadata)
{
    extractBcl(
        bamBlock,
        boost::make_reverse_iterator(randomAccessIt + readMetadata.getLength()),
        &oligo::getReverseBcl,
        readMetadata).base();
    return randomAccessIt + readMetadata.getLength();
}

template <typename RandomAccessIt>
RandomAccessIt extractBcl(
    const BamBlockHeader &bamBlock,
    RandomAccessIt randomAccessIt,
    const flowcell::ReadMetadata &readMetadata)
{
    return bamBlock.isReverse() ?
        extractReverseBcl(bamBlock, randomAccessIt, readMetadata) :
        extractForwardBcl(bamBlock, randomAccessIt, readMetadata);
}


} // namespace bam
} // namespace isaac

#endif // #ifndef iSAAC_BAM_BAM_PARSER_HH
