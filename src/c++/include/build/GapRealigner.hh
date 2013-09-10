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
 ** \file GapRealigner.hh
 **
 ** Attempts to reduce read mismatches by introducing gaps found on other reads.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_GAP_REALIGNER_HH
#define iSAAC_BUILD_GAP_REALIGNER_HH

#include "alignment/Cigar.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "build/BarcodeBamMapping.hh"
#include "build/PackedFragmentBuffer.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "gapRealigner/Gap.hh"
#include "reference/Contig.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace build
{

class RealignerGaps
{
    // All the gaps in the sample
    gapRealigner::Gaps gapGroups_;
    // Deletion gaps sorted by their end position
    gapRealigner::Gaps deletionEndGroups_;

public:
    typedef gapRealigner::Gap GapType;

    void reserve(const size_t gaps);
    void unreserve();

    void addGapsFromFragment(const io::FragmentAccessor &fragment);

    template<typename IteratorT>
    void addGaps(
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
                addGap(gapRealigner::Gap(pos, -decoded.first));
            }
            else if (decoded.second == Cigar::DELETE)
            {
                addGap(gapRealigner::Gap(pos, decoded.first));
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
        const gapRealigner::Gap &gap)
    {
        gapGroups_.push_back(gap);
//        ISAAC_THREAD_CERR_DEV_TRACE("Added " << gaps.back());
    }

    void finalizeGaps();

    size_t getGapsCount() const
    {
        return gapGroups_.size();
    }

    gapRealigner::GapsRange findGaps(
        const unsigned long clusterId,
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition rangeBegin,
        const reference::ReferencePosition rangeEnd,
        gapRealigner::Gaps &foundGaps_) const;
};



/**
 * \brief Attempts to insert gaps found on other fragments while preserving the ones that
 *        are already there.
 */
class GapRealigner
{
    // number of bits that can represent the on/off state for each gap.
    // Currently unsigned short is used to hold the choice
    static const unsigned MAX_GAPS_AT_A_TIME = 10;

    const bool realignGapsVigorously_;
    const bool realignDodgyFragments_;
    const unsigned realignedGapsPerFragment_;
    // Recommended value to be lower than gapOpenCost_ in a way that
    // no less than two mismatches would warrant adding a gap
    const unsigned mismatchCost_;// = 3;
    const unsigned gapOpenCost_;// = 4;
    // Recommended 0 as it does not matter how long the introduced gap is for realignment
    const unsigned gapExtendCost_;// = 0;
    static const int mismatchPercentReductionMin_ = 20;

    const bool clipSemialigned_;

    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics_;
    const std::vector<std::vector<reference::Contig> > &contigList_;

    alignment::Cigar realignedCigars_;
    // In vigorous mode lastAttemptGaps_ is used to avoid realigning against the gaps that have been already tested
    // in the previous pass. TODO: does not seem to improve timing but causes a minor amount of less-well aligned reads
//    gapRealigner::Gaps lastAttemptGaps_;
    gapRealigner::Gaps currentAttemptGaps_;

public:
    typedef gapRealigner::Gap GapType;
    GapRealigner(
        const bool realignGapsVigorously,
        const bool realignDodgyFragments,
        const unsigned realignedGapsPerFragment,
        const unsigned mismatchCost,
        const unsigned gapOpenCost,
        const unsigned gapExtendCost,
        const bool clipSemialigned,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        const std::vector<std::vector<reference::Contig> > &contigList):
            realignGapsVigorously_(realignGapsVigorously),
            realignDodgyFragments_(realignDodgyFragments),
            realignedGapsPerFragment_(realignedGapsPerFragment),
            mismatchCost_(mismatchCost),
            gapOpenCost_(gapOpenCost),
            gapExtendCost_(gapExtendCost),
            clipSemialigned_(clipSemialigned),
            barcodeMetadataList_(barcodeMetadataList),
            barcodeTemplateLengthStatistics_(barcodeTemplateLengthStatistics),
            contigList_(contigList)
    {
//        lastAttemptGaps_.reserve(MAX_GAPS_AT_A_TIME * 10);
        currentAttemptGaps_.reserve(MAX_GAPS_AT_A_TIME * 10);
    }
    GapRealigner(const GapRealigner &that):
        realignGapsVigorously_(that.realignGapsVigorously_),
        realignDodgyFragments_(that.realignDodgyFragments_),
        realignedGapsPerFragment_(that.realignedGapsPerFragment_),
        mismatchCost_(that.mismatchCost_),
        gapOpenCost_(that.gapOpenCost_),
        gapExtendCost_(that.gapExtendCost_),
        clipSemialigned_(that.clipSemialigned_),
        barcodeMetadataList_(that.barcodeMetadataList_),
        barcodeTemplateLengthStatistics_(that.barcodeTemplateLengthStatistics_),
        contigList_(that.contigList_)
    {
//        lastAttemptGaps_.reserve(MAX_GAPS_AT_A_TIME * 10);
        currentAttemptGaps_.reserve(MAX_GAPS_AT_A_TIME * 10);
    }

    void realign(
        RealignerGaps &realignerGaps,
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition binEndPos,
        PackedFragmentBuffer::Index &index,
        io::FragmentAccessor &fragment,
        PackedFragmentBuffer &dataBuffer);

    void reserve(const alignment::BinMetadata& bin)
    {
        // assume each existing cigar gets realignedGapsPerFragment_ gaps introduced...
        realignedCigars_.reserve(bin.getTotalCigarLength() + bin.getTotalElements() * (1 + realignedGapsPerFragment_ * 2));
    }

    void unreserve()
    {
        alignment::Cigar().swap(realignedCigars_);
    }

private:

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
         * \brief   Position of the first base following the last insertion or the first base
         *          that is not part of the last deletion. If there are no indels, equals to beginPos_
         */
        reference::ReferencePosition lastGapEndPos_;
        /*
         * \brief   Position of the base that follows the last non soft-clipped base of the read
         */
        reference::ReferencePosition endPos_;
    };
    friend std::ostream & operator << (std::ostream &os, const GapRealigner::RealignmentBounds &fragmentGaps);

    const gapRealigner::GapsRange findMoreGaps(
        gapRealigner::GapsRange range,
        const gapRealigner::Gaps &gaps,
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition binEndPos);

    const gapRealigner::GapsRange findGaps(
        const unsigned sampleId,
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition binEndPos,
        const reference::ReferencePosition rangeBegin,
        const reference::ReferencePosition rangeEnd);

    bool applyChoice(
        const unsigned short choice,
        const gapRealigner::GapsRange &gaps,
        const reference::ReferencePosition binEndPos,
        const reference::ReferencePosition contigEndPos,
        PackedFragmentBuffer::Index &index,
        const io::FragmentAccessor &fragment);


    struct GapChoice
    {
        GapChoice() : editDistance_(0), mismatches_(0), cost_(0), mappedLength_(0){}
        unsigned editDistance_;
        unsigned mismatches_;
        unsigned cost_;
        unsigned mappedLength_;
    };
    friend std::ostream & operator <<(std::ostream &os, const GapChoice &gapChoice)
    {
        return os << "GapChoice(" << gapChoice.editDistance_ << "ed," << gapChoice.mismatches_ << "mm," << gapChoice.cost_ << "c," << gapChoice.mappedLength_ << "ml)";
    }


    GapChoice verifyGapsChoice(
        const unsigned short choice,
        const gapRealigner::GapsRange &gaps,
        const reference::ReferencePosition newBeginPos,
        const PackedFragmentBuffer::Index &index,
        const io::FragmentAccessor &fragment,
        const std::vector<reference::Contig> &reference);

    bool isBetterChoice(
        const GapChoice &choice,
        const int originalMismatchesPercent,
        const unsigned bestCost,
        const unsigned bestEditDistance) const;

    const RealignmentBounds extractRealignmentBounds(const PackedFragmentBuffer::Index &index);

    bool findStartPos(
        const unsigned short choice,
        const gapRealigner::GapsRange &gaps,
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition binEndPos,
        const PackedFragmentBuffer::Index &index,
        const unsigned pivotGapIndex,
        const reference::ReferencePosition pivotPos,
        reference::ReferencePosition &ret);

    bool compactCigar(
        const std::vector<reference::Contig> &reference,
        const reference::ReferencePosition binEndPos,
        PackedFragmentBuffer::Index &index,
        io::FragmentAccessor &fragment);

    unsigned getAlignmentCost(
        const io::FragmentAccessor &fragment,
        const PackedFragmentBuffer::Index &index,
        unsigned &editDistance,
        int &mismatchesPercent) const;

    void compactRealignedCigarBuffer(
        std::size_t bufferSizeBeforeRealignment,
        PackedFragmentBuffer::Index &index);

    void updatePairDetails(
        const PackedFragmentBuffer::Index &index,
        io::FragmentAccessor &fragment,
        PackedFragmentBuffer &dataBuffer);
};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_GAP_REALIGNER_HH
