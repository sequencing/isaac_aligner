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
 ** \file GapRealigner.hh
 **
 ** Attempts to reduce read mismatches by introducing gaps found on other reads.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_GAP_REALIGNER_HH
#define iSAAC_BUILD_GAP_REALIGNER_HH

#include "alignment/Cigar.hh"
#include "build/BarcodeBamMapping.hh"
#include "build/PackedFragmentBuffer.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "reference/Contig.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace build
{
namespace gapRealigner
{
struct Gap
{
    Gap(const reference::ReferencePosition pos, const int length) :
        pos_(pos), length_(length){}
    /**
     * \brief first position of the indel in reference
     */
    reference::ReferencePosition pos_;
    /**
     * \brief positive value indicates deletion from reference. negative - deletion from data.
     *        zero - just an indication of a position in the reference. Intented to be used for lookups.
     */
    int length_;

    unsigned getLength() const {return std::abs(length_);}
    alignment::Cigar::OpCode getOpCode() const { return isInsertion() ? alignment::Cigar::INSERT : alignment::Cigar::DELETE;}
    bool isInsertion() const {return 0 > length_;}
    bool isDeletion() const {return 0 < length_;}
    bool operator < (const Gap &right) const
    {
/*            const reference::ReferencePosition endPos = (pos_ + length_);
        const reference::ReferencePosition rightEndPos = (right.pos_ + right.length_);
        // out of those that end at the same position, put the longer ones in front as
        // they have a higher chance of overlapping the end of the query range.
        return endPos < rightEndPos ||
            (endPos == rightEndPos && length_ > right.length_);*/

        return pos_ < right.pos_ ||
            (pos_ == right.pos_ && length_ < right.length_);

    }

    bool operator ==(const Gap &right) const
    {
        return pos_ == right.pos_ && length_ == right.length_;
    }

    reference::ReferencePosition getBeginPos() const
    {
        return pos_;
    }
    reference::ReferencePosition getEndPos(const bool fatInserstions) const
    {
        return (isDeletion() || fatInserstions) ? pos_ + std::abs(length_) : pos_;
    }

};

inline std::ostream &operator << (std::ostream &os, const Gap& gap)
{
    return os << "Gap(" << gap.pos_ << "," << gap.length_ << ")";
}
} //namespace gapRealigner


/**
 * \brief Attempts to insert gaps found on other fragments while preserving the ones that
 *        are already there.
 */
class GapRealigner
{
    // number of bits that can represent the on/off state for each gap.
    // Currently unsigned short is used to hold the choice
    static const unsigned MAX_GAPS_AT_A_TIME = 10;

    /// Arbitrarily high number to avoid buffer reallocation
    static const unsigned insertionsPerFragmentMax_ = 1000;

    const bool clipSemialigned_;

    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const std::vector<std::vector<reference::Contig> > &contigList_;
    const BarcodeBamMapping &barcodeBamMapping_;

    typedef std::vector<gapRealigner::Gap> Gaps;

    std::vector<Gaps> sampleGaps_;
    alignment::Cigar realignedCigars_;

public:
    typedef gapRealigner::Gap GapType;
    GapRealigner(
        const bool clipSemialigned,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const std::vector<std::vector<reference::Contig> > &contigList,
        const BarcodeBamMapping &barcodeBamMapping):
            clipSemialigned_(clipSemialigned),
            barcodeMetadataList_(barcodeMetadataList), contigList_(contigList),
            barcodeBamMapping_(barcodeBamMapping),
            sampleGaps_(barcodeBamMapping_.getTotalFiles()){}
    GapRealigner(const GapRealigner &that):
        clipSemialigned_(that.clipSemialigned_),
        barcodeMetadataList_(that.barcodeMetadataList_), contigList_(that.contigList_),
        barcodeBamMapping_(that.barcodeBamMapping_),
        sampleGaps_(that.sampleGaps_){}

    void realign(
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition binEndPos,
        PackedFragmentBuffer::Index &index,
        PackedFragmentBuffer &dataBuffer);

    void realignFast(
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition binEndPos,
        PackedFragmentBuffer::Index &index,
        PackedFragmentBuffer &dataBuffer);

    void reserve(const alignment::BinMetadata& bin);

    void unreserve()
    {
        alignment::Cigar().swap(realignedCigars_);
        BOOST_FOREACH(Gaps& gaps, sampleGaps_)
        {
            Gaps().swap(gaps);
        }
    }

    void addGapsFromFragment(const io::FragmentAccessor &fragment);

    template<typename IteratorT>
    void addGaps(
        const unsigned sampeId,
        reference::ReferencePosition fStrandPosition,
        IteratorT cigarBegin, IteratorT cigarEnd)
    {
        using alignment::Cigar;
        reference::ReferencePosition pos = fStrandPosition;
        IteratorT cigarIterator = cigarBegin;
        for(;cigarEnd != cigarIterator; ++cigarIterator)
        {
            const Cigar::Component decoded = Cigar::decode(*cigarIterator);
            if (decoded.second == Cigar::ALIGN)
            {
                pos += decoded.first;
            }
            else if (decoded.second == Cigar::INSERT)
            {
                addGap(sampeId, gapRealigner::Gap(pos, -decoded.first));
            }
            else if (decoded.second == Cigar::DELETE)
            {
                addGap(sampeId, gapRealigner::Gap(pos, decoded.first));
                pos += decoded.first;
            }
            else if (decoded.second == Cigar::SOFT_CLIP)
            {
                if (fStrandPosition == pos)
                {
                    ISAAC_ASSERT_MSG(cigarBegin == cigarIterator || cigarEnd == cigarIterator + 1,
                                     "First soft clip can be only the first or the last component of the cigar");
                    // not advancing pos as it does not include the soft-clipped area
                }
                else
                {
                    ISAAC_ASSERT_MSG(cigarEnd == cigarIterator + 1,
                                     "At most two soft-clips are expected with second one being the last component of the cigar");
                }
            }
            else
            {
                const boost::format message = boost::format("Unexpected Cigar OpCode: %d") % decoded.second;
                ISAAC_ASSERT_MSG(false, message.str().c_str());
            }
        }
    }

    void addGap(
        const unsigned sampeId,
        const gapRealigner::Gap &gap)
    {
        Gaps &gaps = sampleGaps_.at(sampeId);
        gaps.push_back(gap);
//        ISAAC_THREAD_CERR_DEV_TRACE("Added " << gaps.back());
    }

    void finalizeGaps()
    {
        BOOST_FOREACH(Gaps& gaps, sampleGaps_)
        {
            std::sort(gaps.begin(), gaps.end());
            gaps.erase(std::unique(gaps.begin(), gaps.end()), gaps.end());
        }
    }

    size_t getTotalGapsCount() const
    {
        return std::accumulate(sampleGaps_.begin(), sampleGaps_.end(), 0,
                               boost::bind(std::plus<size_t>(), _1, boost::bind(&std::vector<gapRealigner::Gap>::size, _2)));
    }

private:
    struct GapsRange : public std::pair<Gaps::const_iterator, Gaps::const_iterator>
    {
        typedef std::pair<Gaps::const_iterator, Gaps::const_iterator> BaseType;
        GapsRange(Gaps::const_iterator f, Gaps::const_iterator s): BaseType(f, s) {}
        GapsRange(){}
        bool empty() const {return second == first;}
        unsigned size() const {return second - first;}
    };
    friend std::ostream &operator <<(std::ostream &os, const GapsRange &gaps);

    struct RealignmentBounds
    {
        /*
         * \brief Position of the first non soft-clipped base of the read
         */
        reference::ReferencePosition beginPos_;
        /*
         * \breif   Position of the first insertion base or the first base before the first deletion.
         *          If there are no indels, equals to endPos.
         */
        reference::ReferencePosition firstGapStartPos_;
        /*
         * \brief   Posgiition of the first base following the last insertion or the first base
         *          that is not part of the last deletion. If there are no indels, equals to beginPos_
         */
        reference::ReferencePosition lastGapEndPos_;
        /*
         * \brief   Position of the base that follows the last non soft-clipped base of the read
         */
        reference::ReferencePosition endPos_;
    };
    friend std::ostream & operator << (std::ostream &os, const GapRealigner::RealignmentBounds &fragmentGaps);

    const GapsRange findGaps(
        const unsigned sampleId,
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition rangeBegin,
        const reference::ReferencePosition rangeEnd) const;

    void applyChoice(
        const unsigned short choice,
        const GapsRange &gaps,
        const reference::ReferencePosition binStartPos,
        PackedFragmentBuffer::Index &index,
        const io::FragmentAccessor &fragment,
        PackedFragmentBuffer &dataBuffer);

    unsigned verifyGapsChoice(
        const unsigned short choice,
        const GapsRange &gaps,
        reference::ReferencePosition newBeginPos,
        const reference::ReferencePosition binStartPos,
        const PackedFragmentBuffer::Index &index,
        const io::FragmentAccessor &fragment,
        RealignmentBounds bounds,
        const std::vector<reference::Contig> &reference,
        unsigned &mismatches);
    const RealignmentBounds extractRealignmentBounds(const PackedFragmentBuffer::Index &index);

    bool findStartPos(
        const unsigned short choice,
        const GapsRange &gaps,
        const reference::ReferencePosition binStartPos,
        const PackedFragmentBuffer::Index &index,
        const unsigned pivotGapIndex,
        const reference::ReferencePosition pivotPos,
        reference::ReferencePosition &ret);

    bool compactCigar(
        const reference::ReferencePosition binEndPos,
        PackedFragmentBuffer::Index &index,
        io::FragmentAccessor &fragment);
};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_GAP_REALIGNER_HH
