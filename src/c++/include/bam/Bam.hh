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
 ** \file Bam.hh
 **
 ** \brief collection of helpers for bam serialization.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BAM_BAM_HH
#define iSAAC_BAM_BAM_HH

#include <ostream>

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>

#include "config.h"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace bam
{

struct iTag
{
    iTag(): value_(0)
    {
        tag_[0] = 0;
        tag_[1] = 0;
    }

    iTag(const char tag[2], int value):
        value_(value)
    {
        tag_[0] = tag[0];
        tag_[1] = tag[1];
    }
    char tag_[2];
    static const char val_type_ = 'i';
    int  value_;

    bool empty() const {return !tag_[0];}
    /**
     * \brief returns 0 if tag is not set
     */
    std::size_t size() const { return !empty() ? sizeof(tag_) + sizeof(val_type_) + sizeof(value_) : 0;}
};

struct zTag
{
    zTag() : value_(0), valueEnd_(0)
    {
        tag_[0] = 0;
        tag_[1] = 0;
    }

    zTag(const char tag[2], const char *value):
        value_(value), valueEnd_(value_ ? (value_ + strlen(value_) + 1) : 0)
    {
        tag_[0] = tag[0];
        tag_[1] = tag[1];
    }

    zTag(const char tag[2], const char *value, const char *valueEnd):
        value_(value), valueEnd_(valueEnd)
    {
        tag_[0] = tag[0];
        tag_[1] = tag[1];
    }
    char tag_[2];
    static const char val_type_ = 'Z';
    const char *value_;
    const char *valueEnd_;

    bool empty() const {return !tag_[0];}
    std::size_t size() const {return empty() ? 0 : (sizeof(tag_) + sizeof(val_type_) + std::distance(value_, valueEnd_));}
};

void serialize(std::ostream &os, const char* bytes, std::size_t size);

inline void serialize(std::ostream &os, const char* pStr, const char* pEnd) {
    serialize(os, pStr, std::distance(pStr, pEnd));
}

inline void serialize(std::ostream &os, const char* pStr) {
    serialize(os, pStr, strlen(pStr) + 1);
}

inline void serialize(std::ostream &os, const std::string &str) {
    serialize(os, str.c_str(), str.length() + 1);
}

//todo: provide proper implementation with byte flipping
inline void serialize(std::ostream &os, const int &i) {
    serialize(os, reinterpret_cast<const char*>(&i), sizeof(i));
}

inline void serialize(std::ostream &os, const char &c) {
    serialize(os, &c, sizeof(c));
}

//todo: provide proper implementation with byte flipping
inline void serialize(std::ostream &os, const unsigned &ui) {
    serialize(os, reinterpret_cast<const char*>(&ui), sizeof(ui));
}

inline void serialize(std::ostream &os, const iTag &tag) {
    serialize(os, tag.tag_, sizeof(tag.tag_));
    const char val_type = tag.val_type_;
    serialize(os, val_type);
    serialize(os, tag.value_);
}

inline void serialize(std::ostream &os, const zTag &tag) {
    if (tag.value_)
    {
        serialize(os, tag.tag_, sizeof(tag.tag_));
        const char val_type = tag.val_type_;
        serialize(os, val_type);
        serialize(os, tag.value_, tag.valueEnd_);
    }
}


template <typename T>
void serialize(std::ostream &os, const std::vector<T> &vector) {
    serialize(os, reinterpret_cast<const char*>(&vector.front()), vector.size() * sizeof(T));
}

template <typename IteratorT>
void serialize(std::ostream &os, const std::pair<IteratorT, IteratorT> &pairBeginEnd) {
    serialize(os, reinterpret_cast<const char*>(&*pairBeginEnd.first),
              std::distance(pairBeginEnd.first, pairBeginEnd.second) * sizeof(*pairBeginEnd.first));
}

static const unsigned MAX_LANES_PER_FLOWCELL = 8;
static const unsigned MAX_TILES_PER_LANE = 2048;

template <typename THeader>
void serializeHeader(
    std::ostream &os,
    const std::vector<std::string>& argv,
    const std::string &description,
    const std::vector<std::string>& headerTags,
    const THeader &header)
{
    struct Header
    {
        char magic[4];
        int l_text;
    } __attribute__ ((packed));

    const std::string commandLine(boost::join(argv, " "));

    std::string headerText(
        "@HD\t"
            "VN:1.0\t"
            "SO:coordinate\n"
        "@PG\t"
            "ID:iSAAC\t"
            "PN:iSAAC\t"
            "CL:" + commandLine + "\t" +
            (description.empty() ? std::string() : ("DS:" + description + "\t")) +
            "VN:" + iSAAC_VERSION_FULL +
            "\n");

    BOOST_FOREACH(const std::string &headerTag, headerTags)
    {
        headerText += headerTag + "\n";
    }

    BOOST_FOREACH(const typename THeader::ReadGroupType &readGroup, header.getReadGroups())
    {
        headerText += readGroup.getValue() + "\n";
    }

    const typename THeader::RefSeqsType &refSeqs = header.getRefSequences();

    BOOST_FOREACH(const typename THeader::RefSeqType &refSeq, refSeqs)
    {
        std::string sq = "@SQ\tSN:" + refSeq.name() + "\tLN:" + boost::lexical_cast<std::string>(refSeq.length());

        if (refSeq.bamSqAs().length())
        {
           sq += "\tAS:" + refSeq.bamSqAs();
        }

        if (!refSeq.bamSqUr().empty())
        {
           sq += "\tUR:" + refSeq.bamSqUr();
        }

        if(!refSeq.bamM5().empty())
        {
           sq += "\tM5:" + refSeq.bamM5();
        }

        headerText += sq + "\n";
    }

    Header bamHeader ={ {'B','A','M',1}, int(headerText.size())};

    if (!os.write(reinterpret_cast<char*>(&bamHeader), sizeof(bamHeader))){
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write BAM header into bam stream"));
    }

    // samtools view -H ends up printing the binary zero if it is stored here
    serialize(os, headerText.c_str(), headerText.size());

    serialize(os, boost::numeric_cast<int>(refSeqs.size())); //n_ref
    BOOST_FOREACH(const typename THeader::RefSeqType &refSeq, refSeqs)
    {
        const std::string &name(refSeq.name());
        int l_name(name.length() + 1);
        int l_ref(refSeq.length());
        serialize(os, l_name);
        serialize(os, name);
        serialize(os, l_ref);
    }
}

/*!
  @abstract    Calculate the minimum bin that contains a region [beg,end).
  @param  beg  start of the region, 0-based
  @param  end  end of the region, 0-based
  @return      bin
 */
inline int bam_reg2bin(unsigned beg, unsigned end)
{
    --end;
    if (beg>>14 == end>>14) return 4681 + (beg>>14);
    if (beg>>17 == end>>17) return  585 + (beg>>17);
    if (beg>>20 == end>>20) return   73 + (beg>>20);
    if (beg>>23 == end>>23) return    9 + (beg>>23);
    if (beg>>26 == end>>26) return    1 + (beg>>26);
    return 0;
}

typedef std::vector<flowcell::TileMetadata> TileMetadataList;

template <typename T>
unsigned serializeAlignment(std::ostream &os, T&alignment)
{
    const int refID(alignment.refId());
    const int pos(alignment.pos());

    const char *readName = alignment.readName();
    const std::size_t readNameLength = strlen(readName);
    ISAAC_ASSERT_MSG(0xFF > readNameLength, "Read name length must fit in 8 bit value");

    const unsigned observedLength = alignment.observedLength();
    const unsigned bin_mq_nl(unsigned(bam_reg2bin(pos, pos + (observedLength ? observedLength : 1))) << 16 |
                             unsigned(alignment.mapq()) << 8 |
                             unsigned(readNameLength + 1));

    typedef typename T::CigarBeginEnd CigarBeginEnd;
    const CigarBeginEnd cigarBeginEnd = alignment.cigar();
    const std::size_t cigarLength = std::distance(cigarBeginEnd.first, cigarBeginEnd.second);
    ISAAC_ASSERT_MSG(0xFFFF >= cigarLength, "Cigar length must fit in 16 bit value");

    unsigned flag_nc(unsigned (alignment.flag()) << 16 | (unsigned (cigarLength) & 0xFFFF));
    const int l_seq(alignment.seqLen());
    const int next_RefID(alignment.nextRefId());
    const int next_pos(alignment.nextPos());
    const int tlen(alignment.tlen());

    const std::vector<unsigned char> &seq = alignment.seq();
    const std::vector<unsigned char> &qual = alignment.qual();

    const iTag fragmentSM = alignment.getFragmentSM();
    const iTag fragmentAS = alignment.getFragmentAS();
    const zTag fragmentRG = alignment.getFragmentRG();
    const iTag fragmentNM = alignment.getFragmentNM();
    const zTag fragmentBC = alignment.getFragmentBC();
    const zTag fragmentOC = alignment.getFragmentOC();
    const iTag fragmentZX = alignment.getFragmentZX();
    const iTag fragmentZY = alignment.getFragmentZY();

    const int block_size(  sizeof(refID)
                           + sizeof(pos)
                           + sizeof(bin_mq_nl)
                           + sizeof(flag_nc)
                           + sizeof(l_seq)
                           + sizeof(next_RefID)
                           + sizeof(next_pos)
                           + sizeof(tlen)
                           + readNameLength + 1
                           + cigarLength * sizeof(unsigned)
                           + seq.size()
                           + qual.size()
                           + fragmentSM.size()
                           + fragmentAS.size()
                           + fragmentNM.size()
                           + fragmentBC.size()
                           + fragmentRG.size()
                           + fragmentOC.size()
                           + fragmentZX.size()
                           + fragmentZY.size());

    serialize(os, block_size);
    serialize(os, refID);
    serialize(os, pos);

    serialize(os, bin_mq_nl);
    serialize(os, flag_nc);

    serialize(os, l_seq);
    serialize(os, next_RefID);
    serialize(os, next_pos);
    serialize(os, tlen);

    serialize(os, readName);
    serialize(os, cigarBeginEnd);//4
    serialize(os, seq);  //1
    serialize(os, qual); //2

    if (!fragmentSM.empty()) {serialize(os, fragmentSM);}
    if (!fragmentAS.empty()) {serialize(os, fragmentAS);}

    if (!fragmentRG.empty()) {serialize(os, fragmentRG);}
    if (!fragmentNM.empty()) {serialize(os, fragmentNM);}
    if (!fragmentBC.empty()) serialize(os, fragmentBC);
    if (!fragmentOC.empty()) serialize(os, fragmentOC);

    if (!fragmentZX.empty()) {serialize(os, fragmentZX);}
    if (!fragmentZY.empty()) {serialize(os, fragmentZY);}

    return block_size + sizeof(block_size);
}

void serializeBgzfFooter(std::ostream &os);

} //namespace bam
} // namespace isaac


#endif // iSAAC_BAM_BAM_HH
