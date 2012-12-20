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
 ** \file Fragment.hh
 **
 ** \brief Defines io structures for pre-bam bin fragment io.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_FRAGMENT_HH
#define iSAAC_IO_FRAGMENT_HH

#include "alignment/FragmentMetadata.hh"
#include "alignment/BamTemplate.hh"
#include "alignment/Cigar.hh"

namespace isaac
{
namespace io
{

struct FragmentHeader
{
    FragmentHeader():
        bamTlen_(0),
        observedLength_(0),
        fStrandPosition_(),
        lowClipped_(0),
        highClipped_(0),
        alignmentScore_(0),
        templateAlignmentScore_(0),
        mateFStrandPosition_(),
        readLength_(0),
        cigarLength_(0),
        gapCount_(0),
        editDistance_(0),
        flags_(false, false, false, false, false, false, false, false, false, false),
        tile_(0),
        barcode_(0),
        clusterId_(0)
    {
    }

    FragmentHeader(const alignment::BamTemplate &bamTemplate,
                   const alignment::FragmentMetadata &fragment,
                   const alignment::FragmentMetadata &mate,
                   const unsigned barcodeIdx,
                   const bool mateBinTheSame)
    :
        bamTlen_(getTlen(fragment, mate)),
        observedLength_(fragment.getObservedLength()),
        fStrandPosition_(!fragment.isNoMatch() ?
            fragment.getFStrandReferencePosition() :
            mate.getFStrandReferencePosition()),
        lowClipped_(fragment.lowClipped),
        highClipped_(fragment.highClipped),
        alignmentScore_(fragment.getAlignmentScore()),
        templateAlignmentScore_(
            bamTemplate.isProperPair()
            ? bamTemplate.getAlignmentScore()
            : fragment.getAlignmentScore()),
        mateFStrandPosition_(!mate.isNoMatch() ?
            mate.getFStrandReferencePosition() :
            fragment.getFStrandReferencePosition()),
        readLength_(fragment.getReadLength()),
        cigarLength_(fragment.getCigarLength()),
        gapCount_(fragment.getGapCount()),
        editDistance_(fragment.getEditDistance()),
        flags_(true,
               !fragment.isAligned(),
               !mate.isAligned(),
               fragment.isReverse(),
               mate.isReverse(),
               0 == fragment.getReadIndex(),
               1 == fragment.getReadIndex(),
               !fragment.getCluster().getPf(),
               bamTemplate.isProperPair(),
               mateBinTheSame),
        tile_(fragment.getCluster().getTile()),
        barcode_(barcodeIdx),
        clusterId_(fragment.getCluster().getId())
        //, magic_(magicValue_)
    {
    }

    FragmentHeader(const alignment::BamTemplate &bamTemplate,
                   const alignment::FragmentMetadata &fragment,
                   const unsigned barcodeIdx)
    :
        // According to SAM v1.4 TLEN is 0 for single-ended templates.
        bamTlen_(0),
        observedLength_(fragment.getObservedLength()),
        fStrandPosition_(fragment.getFStrandReferencePosition()),
        lowClipped_(fragment.lowClipped),
        highClipped_(fragment.highClipped),
        alignmentScore_(fragment.getAlignmentScore()),
        templateAlignmentScore_(fragment.getAlignmentScore()),
        mateFStrandPosition_(reference::ReferencePosition::NoMatch),
        readLength_(fragment.getReadLength()),
        cigarLength_(fragment.getCigarLength()),
        gapCount_(fragment.getGapCount()),
        editDistance_(fragment.getEditDistance()),
        flags_(false,
               !fragment.isAligned(),
               true,
               fragment.isReverse(),
               false,
               true,
               true,
               !fragment.getCluster().getPf(),
               false,
               // set mateBinTheSame to true to allow single-ended reads realignment
               true),
        tile_(fragment.getCluster().getTile()),
        barcode_(barcodeIdx),
        clusterId_(fragment.getCluster().getId())
        //, magic_(magicValue_)
    {
    }

    static unsigned getDataLength(const unsigned readLength, const unsigned cigarLength) {
        return readLength + cigarLength * sizeof(alignment::Cigar::value_type);
    }

    unsigned getDataLength() const {
        return getDataLength(readLength_, cigarLength_);
    }

    static unsigned getTotalLength(const unsigned readLength, const unsigned cigarLength) {
        return sizeof(FragmentHeader) + getDataLength(readLength, cigarLength);
    }

    unsigned getTotalLength() const {
        return getTotalLength(readLength_, cigarLength_);
    }

    static unsigned getMaxTotalLength(const unsigned readLength)
    {
        return getTotalLength(readLength, alignment::Cigar::getMaxLength(readLength));
    }

    static unsigned getMinTotalLength(const unsigned readLength)
    {
        return getTotalLength(readLength, alignment::Cigar::getMinLength());
    }

    /*
     * \brief as defined for TLEN in SAM v1.4
     */
    static int getTlen(const reference::ReferencePosition fragmentBeginPos,
                       const reference::ReferencePosition fragmentEndPos,
                       const reference::ReferencePosition mateBeginPos,
                       const reference::ReferencePosition mateEndPos,
                       const bool firstRead)
    {
        // SAM TLEN is by 1 less than the number of bases the template covers
        const unsigned long distance = (
            std::max(fragmentEndPos, mateEndPos).getLocation() -
            std::min(fragmentBeginPos, mateBeginPos).getLocation()) - 1;


        const long ret = fragmentBeginPos < mateBeginPos ? distance :
            (fragmentBeginPos > mateBeginPos || !firstRead) ? -distance : distance ;

//        const long ret = fragmentBeginPos < mateBeginPos ? distance : -distance;

        return ret;
    }

    /*
     * \brief as defined for TLEN in SAM v1.4
     */
    static int getTlen(const alignment::FragmentMetadata &fragment,
                       const alignment::FragmentMetadata &mate)
    {
        return fragment.isAligned() || mate.isAligned() ?
            getTlen(fragment.getBeginReferencePosition(), fragment.getEndReferencePosition(),
                    mate.getBeginReferencePosition(), mate.getEndReferencePosition(), 0 == fragment.getReadIndex()) :
            0;
    }

    bool isAligned() const {return !flags_.unmapped_;}
    /**
     * \brief template length as specified by SAM format
     * TLEN: signed observed Template LENgth. If all segments are mapped to the same reference, the
     * unsigned observed template length equals the number of bases from the leftmost mapped base
     * to the rightmost mapped base. The leftmost segment has a plus sign and the rightmost has a
     * minus sign. The sign of segments in the middle is undefined. It is set as 0 for single-segment
     * template or when the information is unavailable.
     */
    int bamTlen_;


    /**
     * \brief Normal reads      : positive distance between fStrandPosition_ and the base following
     *                            the last unclipped base of the fragment
     *        Shadows           : 0
     *        unaligned cluster : 0
     */
    unsigned observedLength_;


    reference::ReferencePosition fStrandPosition_;

    /// number of bases clipped from begin and end irrespective of alignment
    unsigned short lowClipped_;
    unsigned short highClipped_;

    /**
     * \brief single fragment alignment score
     */
    unsigned short alignmentScore_;

    /**
     * \brief template alignment score:
     *          for normal pairs:     pair alignment score.
     *          for singleton/shadow: singleton alignment score
     *          for others:           0
     */
    unsigned short templateAlignmentScore_;

    /**
     * \brief forward-strand position of mate
     */
    reference::ReferencePosition mateFStrandPosition_;

    /**
     * \brief number of nucleotides in the fragment
     */
    unsigned short readLength_;

    /**
     * \brief number of operations in the CIGAR
     */
    unsigned short cigarLength_;

    /**
     * \brief number insertion and deletion operations in the CIGAR
     */
    unsigned short gapCount_;
    /**
     * \brief edit distance
     */
    unsigned short editDistance_;

    struct Flags
    {
        Flags(bool paired, bool unmapped, bool mateUnmapped, bool reverse, bool mateReverse,
              bool firstRead, bool secondRead, bool failFilter, bool properPair, bool mateBinTheSame):
            paired_(paired), unmapped_(unmapped), mateUnmapped_(mateUnmapped),
            reverse_(reverse), mateReverse_(mateReverse), firstRead_(firstRead), secondRead_(secondRead),
            failFilter_(failFilter), properPair_(properPair), mateBinTheSame_(mateBinTheSame), duplicate_(false){}
        bool paired_ : 1;
        bool unmapped_ : 1;
        bool mateUnmapped_ : 1;
        bool reverse_ : 1;
        bool mateReverse_ : 1;
        bool firstRead_ : 1;
        bool secondRead_ : 1;
        bool failFilter_ : 1;
        bool properPair_ : 1;
        /// true if mate is located in the same bin. Allows for template-changing operations (such as realignement) during bin-sort
        bool mateBinTheSame_ : 1;
        bool duplicate_ : 1;
    } flags_;

    /**
     * \brief 0-based unique tile index
     */
    unsigned long tile_;

    /**
     * \brief 0-based unique barcode index
     */
    unsigned long barcode_;

    /**
     * \brief 0-based cluster index in the tile
     */
    unsigned long clusterId_;

//    unsigned short magic_;
//    static const unsigned short magicValue_ = 0xb1a;

    /// \return number of bases clipped on the left side of the fragment in respect to the reference
    unsigned short leftClipped() const {return flags_.reverse_ ? highClipped_ : lowClipped_;}
    /// \return number of bases clipped on the right side of the fragment in respect to the reference
    unsigned short rightClipped() const {return flags_.reverse_ ? lowClipped_ : highClipped_;}
};

inline std::ostream & operator <<(std::ostream &os, const FragmentHeader::Flags &headerFlags)
{
    return os << "FragmentHeader::Flags(" <<
        (headerFlags.paired_ ? "pe" : "") << "|" <<
        (headerFlags.unmapped_ ? "u":"") << "|" <<
        (headerFlags.mateUnmapped_ ? "mu" : "") << "|" <<
        (headerFlags.reverse_ ? "r" : "") << "|" <<
        (headerFlags.mateReverse_ ? "mr" : "") << "|" <<
        (headerFlags.firstRead_ ? "r1" : "") << "|" <<
        (headerFlags.secondRead_ ? "r2" : "") << "|" <<
        (headerFlags.failFilter_ ? "ff" : "") << "|" <<
        headerFlags.properPair_ << "|" <<
        headerFlags.mateBinTheSame_ << //"," <<
        ")";
}

inline std::ostream & operator <<(std::ostream &os, const FragmentHeader &header)
{
    return os << "FragmentHeader(" <<
        header.bamTlen_ << "," <<
        header.fStrandPosition_ << "," <<
        header.lowClipped_ << ":" << header.highClipped_ << "lchc," <<
        header.alignmentScore_ << "," <<
        header.templateAlignmentScore_ << "," <<
        header.mateFStrandPosition_ << "," <<
        header.readLength_ << "rl," <<
        header.observedLength_ << "ol," <<
        header.cigarLength_ << "cl," <<
        header.gapCount_ << "g," <<
        header.editDistance_ << "ed," <<
        header.flags_ << "," <<
        header.tile_ << "," <<
        header.barcode_ << "," <<
        header.clusterId_ << "id," <<
        header.leftClipped() << "lc," <<
        header.rightClipped() << "rc," <<
        ")";
}


struct FragmentAccessor : public FragmentHeader
{
    const unsigned char *basesBegin() const {
        return reinterpret_cast<const unsigned char*>(this) + sizeof(FragmentHeader);
    }

    const unsigned char *unmaskedBasesBegin() const
    {
        if (cigarLength_)
        {
            using alignment::Cigar;
            const Cigar::Component decoded = Cigar::decode(*cigarBegin());
            if (decoded.first == Cigar::SOFT_CLIP)
            {
                return basesBegin() + decoded.second;
            }
        }
        return basesBegin();
    }
    const unsigned char *basesEnd() const { return basesBegin() + readLength_; }

    const unsigned *cigarBegin() const { return reinterpret_cast<const unsigned*>(basesEnd()); }
    const unsigned *cigarEnd() const { return cigarBegin() + cigarLength_; }
protected:
    // this structure cannot be directly instantiated. Use it to cast the byte pointers to access
    // data read or mapped from a binary file
    FragmentAccessor(){}
};

inline std::ostream & operator <<(std::ostream &os, const FragmentAccessor &fragment)
{
    const FragmentHeader &header = fragment;
    return os << "FragmentAccessor(" << header <<
        ", " << alignment::Cigar::toString(fragment.cigarBegin(), fragment.cigarEnd()) <<
        ")";
}

} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_FRAGMENT_INDEX_HH
