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
 ** \file GapRealigner.cpp
 **
 ** Attempts to reduce read mismatches by introducing gaps found on other reads.
 **
 ** \author Roman Petrovski
 **/
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "alignment/BandedSmithWaterman.hh"
#include "build/GapRealigner.hh"

#include "SemialignedEndsClipper.hh"

namespace isaac
{
namespace build
{

inline std::ostream & operator << (std::ostream &os, const GapRealigner::RealignmentBounds &fragmentGaps)
{
    return os << "FragmentGaps(" <<
        fragmentGaps.beginPos_ << "," <<
        fragmentGaps.firstGapStartPos_ << "," <<
        fragmentGaps.lastGapEndPos_ << "," <<
        fragmentGaps.endPos_ << ")";
}

void GapRealigner::reserve(const alignment::BinMetadata& bin)
{
    std::vector<size_t> gapsBySample(barcodeBamMapping_.getTotalFiles(), 0);
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList_)
    {
        gapsBySample.at(barcodeBamMapping_.getFileIndex(barcode)) +=
            bin.getBarcodeGapCount(barcode.getIndex());
    }
    unsigned sampleId = 0;
    BOOST_FOREACH(const size_t sampleGaps, gapsBySample)
    {
        sampleGaps_.at(sampleId++).reserve(sampleGaps);
    }
    // assume each cigar is one operation long and gets one indel introduced in the middle...
    realignedCigars_.reserve(bin.getTotalCigarLength() * 3);
}

void GapRealigner::addGapsFromFragment(const io::FragmentAccessor &fragment)
{
    if (fragment.gapCount_)
    {
        const unsigned sampleId = barcodeBamMapping_.getFileIndex(fragment.barcode_);
        addGaps(sampleId, fragment.fStrandPosition_, fragment.cigarBegin(), fragment.cigarEnd());
    }
}

/*
inline std::ostream &operator <<(std::ostream &os, const std::pair<Gaps::const_iterator, Gaps::const_iterator> &gaps)
{
    if (gaps.second == gaps.first)
    {
        return os << "(no gaps)";
    }

    BOOST_FOREACH(const gapRealigner::Gap &gap, gaps)
    {
        os << gap << ",";
    }
    return os;
}
*/

inline std::ostream &operator <<(std::ostream &os, const GapRealigner::GapsRange &gaps)
{
    if (gaps.second == gaps.first)
    {
        return os << "(no gaps)";
    }

    BOOST_FOREACH(const gapRealigner::Gap &gap, std::make_pair(gaps.first, gaps.second))
    {
        os << gap << ",";
    }
    return os;
}

/**
 * \brief Find gaps given the position range.
 */
const GapRealigner::GapsRange GapRealigner::findGaps(
    const unsigned sampleId,
    const reference::ReferencePosition binStartPos,
    const reference::ReferencePosition rangeBegin,
    const reference::ReferencePosition rangeEnd) const
{
    GapsRange ret;
    const reference::ReferencePosition earliestPossibleGapBegin =
        rangeBegin < (binStartPos + alignment::BandedSmithWaterman::widestGapSize)  ?
            binStartPos : rangeBegin - alignment::BandedSmithWaterman::widestGapSize;

    //ISAAC_THREAD_CERR << "findGaps effective range: [" << earliestPossibleGapBegin << ";" << rangeEnd << ")" << std::endl;

    const Gaps &gaps = sampleGaps_.at(sampleId);
    ISAAC_THREAD_CERR_DEV_TRACE("findGaps all gaps: " << GapRealigner::GapsRange(gaps.begin(), gaps.end()));
    // the first one that ends on or after the rangeBegin
    ret.first = std::lower_bound(gaps.begin(),
                                 gaps.end(),
                                 gapRealigner::Gap(earliestPossibleGapBegin, -1000000));
    // the first one that ends on or after the rangeEnd
    ret.second = std::lower_bound(ret.first, gaps.end(), gapRealigner::Gap(rangeEnd, 0));
    return ret;
}

const GapRealigner::RealignmentBounds GapRealigner::extractRealignmentBounds(
    const PackedFragmentBuffer::Index &index)
{
    using alignment::Cigar;
    RealignmentBounds ret = {index.pos_, index.pos_, index.pos_, index.pos_};

    const unsigned* cigarIterator = index.cigarBegin_;
    for(;index.cigarEnd_ != cigarIterator; ++cigarIterator)
    {
        const Cigar::Component decoded = Cigar::decode(*cigarIterator);
        if (decoded.second == Cigar::ALIGN)
        {
            ret.firstGapStartPos_ += decoded.first;
            ret.endPos_ += decoded.first;
        }
        else if (decoded.second == Cigar::INSERT)
        {
            ret.lastGapEndPos_ = ret.endPos_;
            ++cigarIterator;
            break;
        }
        else if (decoded.second == Cigar::DELETE)
        {
            ret.endPos_ += decoded.first;
            ret.lastGapEndPos_ = ret.endPos_;
            ++cigarIterator;
            break;
        }
        else if (decoded.second == Cigar::SOFT_CLIP)
        {
            if (index.cigarBegin_ == cigarIterator)
            {
                ISAAC_ASSERT_MSG(ret.firstGapStartPos_ == ret.beginPos_, "first soft clip must happen before anything");
                ISAAC_ASSERT_MSG(ret.lastGapEndPos_ == ret.beginPos_, "first soft clip must happen before anything");
                ISAAC_ASSERT_MSG(ret.endPos_ == ret.beginPos_, "first soft clip must happen before anything");
            }
            else
            {
                ISAAC_ASSERT_MSG(index.cigarEnd_ == cigarIterator + 1,
                                 "At most two soft-clips are expected with second one being the last component of the cigar");
            }
        }
        else
        {
            const boost::format message = boost::format("Unexpected Cigar OpCode: %d") % decoded.second;
            ISAAC_ASSERT_MSG(false, message.str().c_str());
        }
    }

    for(;index.cigarEnd_ != cigarIterator; ++cigarIterator)
    {
        const Cigar::Component decoded = Cigar::decode(*cigarIterator);
        if (decoded.second == Cigar::ALIGN)
        {
            ret.endPos_ += decoded.first;
        }
        else if (decoded.second == Cigar::INSERT)
        {
            ret.lastGapEndPos_ = ret.endPos_;
        }
        else if (decoded.second == Cigar::DELETE)
        {
            ret.endPos_ += decoded.first;
            ret.lastGapEndPos_ = ret.endPos_;
        }
        else if (decoded.second == Cigar::SOFT_CLIP)
        {
            ISAAC_ASSERT_MSG(index.cigarEnd_ == cigarIterator + 1,
                             "Last soft-clip has to be the last component of the cigar");
        }
        else
        {
            const boost::format message = boost::format("Unexpected Cigar OpCode: %d") % decoded.second;
            ISAAC_ASSERT_MSG(false, message.str().c_str());
        }
    }

    ISAAC_ASSERT_MSG(ret.endPos_ >= ret.lastGapEndPos_, "Fragment gap cannot end outside the fragment");
    return ret;
}

static unsigned countMismatches(
    const std::vector<reference::Contig> &reference,
    const unsigned char *basesIterator,
    const reference::ReferencePosition pos,
    unsigned length)
{
    const reference::Contig &contig = reference.at(pos.getContigId());
    std::vector<char>::const_iterator referenceBaseIt = contig.forward_.begin() + pos.getPosition();
    unsigned mismatches = 0;
    BOOST_FOREACH(const unsigned char readBase, std::make_pair(basesIterator, basesIterator + length))
    {
        if (contig.forward_.end() == referenceBaseIt)
        {
            break;
        }
        mismatches += *referenceBaseIt++ != oligo::getBase(0x3 & readBase);
    }

    referenceBaseIt = contig.forward_.begin() + pos.getPosition();
/*    ISAAC_THREAD_CERR << mismatches << " mismatches read '" <<
        oligo::bclToString(basesIterator, length) << "' ref '" <<
        std::string(referenceBaseIt,
                    referenceBaseIt +
                    std::min<unsigned>(length, std::distance(referenceBaseIt, contig.forward_.end()))) << "'" <<
                    std::endl;*/

    return mismatches;
}

/**
 * \brief Adjusts template length and mate fields
 *
 * \param beginPosShift change in the fragment alignment position. negative for insertions, positive for deletions
 * \param endPosShift   change in the fragment end alignment position. negative for insertions, positive for deletions
 */
static void updatePairDetails(
    const PackedFragmentBuffer::Index &index,
    io::FragmentAccessor &fragment,
    PackedFragmentBuffer &dataBuffer)
{
    if (!index.hasMate() || fragment.flags_.mateUnmapped_)
    {
        const int oldBamTlen = fragment.bamTlen_;
        fragment.bamTlen_ = fragment.bamTlen_ < 0 ?
            (-fragment.observedLength_ + 1) : (fragment.observedLength_ - 1);
        if (index.hasMate())
        {
            io::FragmentAccessor &mate  = dataBuffer.getMate(index);
            ISAAC_ASSERT_MSG2(mate.bamTlen_ == -oldBamTlen,
                              "Incoherent BAM template length between fragment and unmapped mate %s %s mate %s", index % fragment % mate);
            mate.bamTlen_ = -fragment.bamTlen_;
            fragment.mateFStrandPosition_ = fragment.fStrandPosition_;
            mate.mateFStrandPosition_ = fragment.fStrandPosition_;
            mate.fStrandPosition_ = fragment.fStrandPosition_;
        }
        return;
    }

    io::FragmentAccessor &mate  = dataBuffer.getMate(index);

    ISAAC_ASSERT_MSG2(mate.bamTlen_ == -fragment.bamTlen_,
                      "Incoherent BAM template length between fragment and mate %s %s mate %s", index % fragment % mate);

    const reference::ReferencePosition fragmentBeginPos = fragment.fStrandPosition_;
    const reference::ReferencePosition fragmentEndPos = fragmentBeginPos + fragment.observedLength_;
    const reference::ReferencePosition mateBeginPos = fragment.mateFStrandPosition_;
    const reference::ReferencePosition mateEndPos = mateBeginPos + mate.observedLength_;

    fragment.bamTlen_ = io::FragmentAccessor::getTlen(fragmentBeginPos, fragmentEndPos, mateBeginPos, mateEndPos, fragment.flags_.firstRead_);
    mate.bamTlen_ = -fragment.bamTlen_;
    mate.mateFStrandPosition_ = fragment.fStrandPosition_;
}

/**
 * \breif Collapses gaps on the ends of CIGAR into soft-clips
 *        This keeps CASAVA happy and probably in general is a right thing to do.
 *        Also adjusts fragment position and observed length.
 *
 * \return false, when the compacting the cigar would move the read past the binEndPos.
 *         If false is returned index and fragment are guaranteed to be unchanged.
 */
bool GapRealigner::compactCigar(
    const reference::ReferencePosition binEndPos,
    PackedFragmentBuffer::Index &index,
    io::FragmentAccessor &fragment)
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, " Compacting " << fragment);

    ISAAC_ASSERT_MSG2(alignment::Cigar::getReadLength(index.cigarBegin_, index.cigarEnd_) == fragment.readLength_,
                      "Broken CIGAR before compacting: %s %s %s",
                      alignment::Cigar::toString(index.cigarBegin_, index.cigarEnd_) %
                      index % fragment);

    unsigned short newEditDistance = fragment.editDistance_;
    using alignment::Cigar;
    const uint32_t *cigarIterator = index.cigarBegin_;
    unsigned softClipStart = 0;
    bool needCompacting = false;
    for (; index.cigarEnd_ != cigarIterator; ++cigarIterator)
    {
        Cigar::Component decoded = Cigar::decode(*cigarIterator);
        if (Cigar::ALIGN == decoded.second)
        {
            break;
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
            if (index.cigarBegin_ != cigarIterator)
            {
                ISAAC_ASSERT_MSG(index.cigarEnd_ == cigarIterator + 1,
                                 "At most two soft-clips are expected with second one being the last component of the cigar");
            }
            softClipStart += decoded.first;
        }
        else if (Cigar::INSERT == decoded.second)
        {
            needCompacting = true;
            softClipStart += decoded.first;
            newEditDistance -= decoded.first;
        }
        else if (Cigar::DELETE == decoded.second)
        {
            needCompacting = true;
            if (binEndPos <= index.pos_ + decoded.first)
            {
                return false;
            }
            index.pos_ += decoded.first;
            newEditDistance -= decoded.first;
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation");
        }
    }

    const uint32_t *cigarBackwardsIterator = index.cigarEnd_ - 1;
    unsigned softClipEnd = 0;
    for (;cigarIterator != cigarBackwardsIterator; --cigarBackwardsIterator)
    {
        Cigar::Component decoded = Cigar::decode(*cigarBackwardsIterator);
        if (Cigar::ALIGN == decoded.second)
        {
            break;
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
            ISAAC_ASSERT_MSG(index.cigarEnd_ == cigarBackwardsIterator + 1,
                             "At most two soft-clips are expected with second one being the last component of the cigar");
            softClipEnd += decoded.first;
        }
        else if (Cigar::INSERT == decoded.second)
        {
            needCompacting = true;
            softClipEnd += decoded.first;
            newEditDistance -= decoded.first;
        }
        else if (Cigar::DELETE == decoded.second)
        {
            needCompacting = true;
            newEditDistance -= decoded.first;
            // nothing to be done for trailing deletes;
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation");
        }
    }

    if (needCompacting)
    {
        size_t before = realignedCigars_.size();

        if (softClipStart)
        {
            ISAAC_ASSERT_MSG(realignedCigars_.capacity() > realignedCigars_.size(), "Realigned CIGAR buffer is out of capacity");
            realignedCigars_.push_back(Cigar::encode(softClipStart, Cigar::SOFT_CLIP));
        }
        ISAAC_ASSERT_MSG2 (std::distance(cigarIterator, cigarBackwardsIterator + 1) < 200 &&
                           std::distance(cigarIterator, cigarBackwardsIterator + 1) > 0,
                           "Suspiciously long %d CIGAR in the middle of compacting: %s %s %s",
                           std::distance(cigarIterator, cigarBackwardsIterator + 1) %
                           alignment::Cigar::toString(index.cigarBegin_, index.cigarEnd_) %
                           index % fragment);

        ISAAC_ASSERT_MSG2(alignment::Cigar::getReadLength(cigarIterator, cigarBackwardsIterator + 1) + softClipStart + softClipEnd == fragment.readLength_,
                          "Broken CIGAR in the middle of compacting: %s %s %s",
                          alignment::Cigar::toString(cigarIterator, cigarBackwardsIterator + 1) %
                          index % fragment);
        ISAAC_ASSERT_MSG(realignedCigars_.capacity() >= (realignedCigars_.size() + std::distance(cigarIterator, cigarBackwardsIterator + 1)), "Realigned CIGAR buffer is out of capacity");
        realignedCigars_.insert(realignedCigars_.end(), cigarIterator, cigarBackwardsIterator + 1);
        if (softClipEnd)
        {
            ISAAC_ASSERT_MSG(realignedCigars_.capacity() > realignedCigars_.size(), "Realigned CIGAR buffer is out of capacity");
            realignedCigars_.push_back(Cigar::encode(softClipEnd, Cigar::SOFT_CLIP));
        }

        index.cigarBegin_ = &realignedCigars_.at(before);
        index.cigarEnd_ = &realignedCigars_.back() + 1;
    }

    reference::ReferencePosition newEndPos = index.pos_;
    for (;cigarIterator != cigarBackwardsIterator + 1; ++cigarIterator)
    {
        Cigar::Component decoded = Cigar::decode(*cigarIterator);
        if (Cigar::ALIGN == decoded.second)
        {
            newEndPos += decoded.first;
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
            ISAAC_ASSERT_MSG(false, "At most two soft-clips are expected. Both at the ends of the CIGAR");
        }
        else if (Cigar::INSERT == decoded.second)
        {
            // do nothing. they don't affect positions
        }
        else if (Cigar::DELETE == decoded.second)
        {
            newEndPos += decoded.first;
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation");
        }
    }

    ISAAC_ASSERT_MSG2(alignment::Cigar::getReadLength(index.cigarBegin_, index.cigarEnd_) == fragment.readLength_,
                      "Broken CIGAR after compacting: %s %s %s",
                      alignment::Cigar::toString(index.cigarBegin_, index.cigarEnd_) %
                      index % fragment);

    fragment.editDistance_ = newEditDistance;
    fragment.fStrandPosition_ = index.pos_;
    fragment.observedLength_ = newEndPos - fragment.fStrandPosition_;

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, " Compacted " << fragment);

    return true;
}

/**
 * \brief bits in choice determine whether the corresponding gaps are on or off
 *
 * \return mismatches or -1U if choice is inapplicable.
 */
unsigned GapRealigner::verifyGapsChoice(
    const unsigned short choice,
    const GapsRange &gaps,
    reference::ReferencePosition newBeginPos,
    const reference::ReferencePosition binStartPos,
    const PackedFragmentBuffer::Index &index,
    const io::FragmentAccessor &fragment,
    RealignmentBounds bounds,
    const std::vector<reference::Contig> &reference,
    unsigned &editDistance)
{
    using alignment::Cigar;

    int basesLeft = fragment.readLength_;
//    ISAAC_THREAD_CERR << "GapRealigner::verifyGapsChoice basesLeft=" << basesLeft << std::endl;
/*
    Cigar::Component decoded = Cigar::decode(*index.cigarBegin_);
    if (Cigar::SOFT_CLIP == decoded.second)
    {
        if (newBeginPos - binStartPos > decoded.first)
        {
            ISAAC_THREAD_CERR << "tada" << std::endl;

            newBeginPos -= decoded.first;
        }
        else
        {
            ISAAC_THREAD_CERR << "tudu" << std::endl;
            basesLeft -= decoded.first;
        }
    }
*/

    reference::ReferencePosition lastGapEndPos = newBeginPos;
    reference::ReferencePosition lastGapBeginPos;
    bool lastGapWasInsertion = false;
    unsigned currentGapIndex = 0;
    unsigned mismatches = 0;
    BOOST_FOREACH(const gapRealigner::Gap& gap, std::make_pair(gaps.first, gaps.second))
    {
        if (choice & (1 << currentGapIndex))
        {
            if (gap.getEndPos(true) <= lastGapEndPos)// || gap.getBeginPos() > lastGapEndPos + basesLeft)
            {
                // the choice requires a gap that cannot be applied.
                // just bail out. there will be another choice just like
                // this one but without the useless gap
//                ISAAC_THREAD_CERR << gap << "does not fit in range (" << lastGapEndPos <<
//                    ";" << lastGapEndPos + basesLeft << ")" << std::endl;
                return -1U;
            }
//            ISAAC_THREAD_CERR << " lastGapEndPos=" << lastGapEndPos << std::endl;

            const reference::ReferencePosition gapClippedBeginPos = std::max(gap.getBeginPos(), newBeginPos);

            if (gapClippedBeginPos < lastGapEndPos)
            {
                return -1U;
                // Allowing overlapping deletions is tricky because it is hard to track back the
                // newBeginPos from the pivot see findStartPos.
            }


            if (gapClippedBeginPos == lastGapBeginPos && (!lastGapWasInsertion || !gap.isInsertion()))
            {
                // The only case where it makes sense to allow two or more gaps starting at the same
                // position is when we want to combine multiple insertions into a larger
                // one. Other cases:
                // deletion/deletion - is disallowed above
                // insertion/deletion - does not make sense (and cause trouble SAAC-253)
                return -1U;
            }

            if (gap.getBeginPos() < newBeginPos + fragment.leftClipped())
            {
                // gap begins inside alignment-independent clipping on the left. Ignore this choice of gaps
                return -1U;
            }


//            ISAAC_THREAD_CERR << " gapClippedBeginPos=" << gapClippedBeginPos << std::endl;

            const long mappedBases = std::min<unsigned>(basesLeft, gapClippedBeginPos - lastGapEndPos);
//            ISAAC_THREAD_CERR << " mappedBases=" << mappedBases << " basesLeft=" << basesLeft << std::endl;

            const unsigned mm = countMismatches(reference,
                                                fragment.basesBegin() + (fragment.readLength_ - basesLeft),
                                                lastGapEndPos, mappedBases);
            editDistance += mm;
            mismatches += mm;
            basesLeft -= mappedBases;

            if (gap.isInsertion())
            {
                const unsigned clippedGapLength =
                    std::min<unsigned>(basesLeft, gap.getEndPos(true) - gapClippedBeginPos);
//                ISAAC_THREAD_CERR << " clippedGapLength=" << clippedGapLength << std::endl;
                editDistance += clippedGapLength;
                basesLeft -= clippedGapLength;
                // insertions referenece end pos is same as begin pos
                lastGapEndPos = std::max(gapClippedBeginPos, gap.getBeginPos());
                lastGapWasInsertion = true;
            }
            else
            {
                const unsigned clippedGapLength = gap.getEndPos(true) - gapClippedBeginPos;
                //ISAAC_THREAD_CERR << "clippedGapLength=" << clippedGapLength << std::endl;
                editDistance += clippedGapLength;
                lastGapEndPos = gap.getEndPos(false);
                lastGapWasInsertion = false;
            }
            lastGapBeginPos = gapClippedBeginPos;

            if (basesLeft < fragment.rightClipped())
            {
                // Gap overlaps alignment-independent clipping on the right. Ignore this choice of gaps
                return -1U;
            }

        }
        ++currentGapIndex;
    }
    if(basesLeft)
    {
        ISAAC_ASSERT_MSG2(0 < basesLeft, "Spent more than readLength. basesLeft: %d choice: %d %s %s, "
            "newBeginPos %s lastGapEndPos %s gaps: %s",
            basesLeft % int(choice) % index % fragment %
            newBeginPos % lastGapEndPos % gaps);
        const unsigned mm = countMismatches(reference, fragment.basesBegin() + (fragment.readLength_ - basesLeft), lastGapEndPos, basesLeft);
        editDistance += mm;
        mismatches += mm;
    }

    return mismatches;

}

void GapRealigner::applyChoice(
    const unsigned short choice,
    const GapsRange &gaps,
    const reference::ReferencePosition binStartPos,
    PackedFragmentBuffer::Index &index,
    const io::FragmentAccessor &fragment,
    PackedFragmentBuffer &dataBuffer)
{
    reference::ReferencePosition newBeginPos = index.pos_;
    using alignment::Cigar;

    const size_t before = realignedCigars_.size();

    int basesLeft = fragment.readLength_;
//    ISAAC_THREAD_CERR << "GapRealigner::applyChoice basesLeft=" << basesLeft << std::endl;

/*
    Cigar::Component decoded = Cigar::decode(*index.cigarBegin_);
    if (Cigar::SOFT_CLIP == decoded.second)
    {
        ISAAC_THREAD_CERR << "Soft clip in the original CIGAR newBeginPos=" << newBeginPos  << "binStartPos=" << binStartPos << std::endl;
        if (newBeginPos - binStartPos > decoded.first)
        {
            newBeginPos -= decoded.first;
        }
        else
        {
            basesLeft -= decoded.first;
            ISAAC_ASSERT_MSG(realignedCigars_.capacity() > realignedCigars_.size(), "Realigned CIGAR buffer is out of capacity");
            realignedCigars_.push_back(Cigar::encode(decoded.first, Cigar::SOFT_CLIP));
        }
        ISAAC_THREAD_CERR << "Soft clip in the original CIGAR done newBeginPos=" << newBeginPos  << "binStartPos=" << binStartPos << std::endl;
    }
*/

    if (fragment.leftClipped())
    {
        Cigar::Component decoded = Cigar::decode(*fragment.cigarBegin());
        ISAAC_ASSERT_MSG(Cigar::SOFT_CLIP == decoded.second, "Original CIGAR is expected to have soft clip at the start");
        ISAAC_ASSERT_MSG(fragment.leftClipped() <= decoded.first, "Original CIGAR soft clip at the start is shorter than left-clipped bases");
        realignedCigars_.push_back(Cigar::encode(fragment.leftClipped(), Cigar::SOFT_CLIP));
        newBeginPos += fragment.leftClipped();
        basesLeft -= fragment.leftClipped();
    }

    if (fragment.rightClipped())
    {
        Cigar::Component decoded = Cigar::decode(*(fragment.cigarEnd() - 1));
        ISAAC_ASSERT_MSG(Cigar::SOFT_CLIP == decoded.second, "Original CIGAR is expected to have soft clip at the end");
        ISAAC_ASSERT_MSG(fragment.rightClipped() <= decoded.first, "Original CIGAR soft clip at the end is shorter than right-clipped bases");
        basesLeft -= fragment.rightClipped();
        // actual cigar is appended at the end of the function
    }

    reference::ReferencePosition lastGapEndPos = newBeginPos;
    unsigned currentGapIndex = 0;
    Cigar::OpCode lastOperation = Cigar::UNKNOWN;
    BOOST_FOREACH(const gapRealigner::Gap& gap, std::make_pair(gaps.first, gaps.second))
    {
        if (choice & (1 << currentGapIndex))
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, " Applying " << gap << " to " << fragment);

//            ISAAC_THREAD_CERR << "lastGapEndPos=" << lastGapEndPos << std::endl;

            const reference::ReferencePosition gapClippedBeginPos = std::max(gap.getBeginPos(), newBeginPos);
//            ISAAC_THREAD_CERR << "gapClippedBeginPos=" << gapClippedBeginPos << std::endl;

            if (gapClippedBeginPos >= lastGapEndPos)
            {
                const long mappedBases = std::min<unsigned>(basesLeft, gapClippedBeginPos - lastGapEndPos);
                if (mappedBases)
                {
                    // avoid 0M in CIGAR
                    ISAAC_ASSERT_MSG(realignedCigars_.capacity() > realignedCigars_.size(), "Realigned CIGAR buffer is out of capacity");
                    realignedCigars_.push_back(Cigar::encode(mappedBases, Cigar::ALIGN));
                    basesLeft -= mappedBases;
                }
//                ISAAC_THREAD_CERR << "mappedBases=" << mappedBases << std::endl;

                if (gap.isInsertion())
                {
                    const unsigned clippedGapLength =
                        std::min<unsigned>(basesLeft, gap.getEndPos(true) - gapClippedBeginPos);
                    //ISAAC_THREAD_CERR << "clippedGapLength=" << clippedGapLength << std::endl;
                    if (alignment::Cigar::INSERT == lastOperation && !mappedBases)
                    {
                        //GATK does not like 2I1I type cigars
                        realignedCigars_.back() =
                            Cigar::encode(Cigar::decode(realignedCigars_.back()).first + clippedGapLength,
                                          alignment::Cigar::INSERT);
                    }
                    else
                    {
                        realignedCigars_.push_back(Cigar::encode(clippedGapLength, gap.getOpCode()));
                        ISAAC_ASSERT_MSG(realignedCigars_.capacity() > realignedCigars_.size(), "Realigned CIGAR buffer is out of capacity");
                        lastOperation = alignment::Cigar::INSERT;
                    }

                    basesLeft -= clippedGapLength;
                    lastGapEndPos = std::max(gapClippedBeginPos, gap.getBeginPos());
//                    ISAAC_THREAD_CERR << "2nd lastGapEndPos=" << lastGapEndPos << std::endl;
                }
                else
                {
                    const unsigned clippedGapLength = gap.getEndPos(true) - gapClippedBeginPos;
//                    ISAAC_THREAD_CERR << "clippedGapLength=" << clippedGapLength << std::endl;
                    if (alignment::Cigar::DELETE == lastOperation && !mappedBases)
                    {
                        //assuming GATK does not like 2D1D type cigars either...
                        realignedCigars_.back() =
                            Cigar::encode(Cigar::decode(realignedCigars_.back()).first + clippedGapLength,
                                          alignment::Cigar::DELETE);
                    }
                    else
                    {
                        ISAAC_ASSERT_MSG(realignedCigars_.capacity() > realignedCigars_.size(), "Realigned CIGAR buffer is out of capacity");
                        realignedCigars_.push_back(Cigar::encode(clippedGapLength, gap.getOpCode()));
                        lastOperation = alignment::Cigar::DELETE;
                    }

                    lastGapEndPos = gap.getEndPos(false);
//                    ISAAC_THREAD_CERR << "2nd lastGapEndPos=" << lastGapEndPos << std::endl;
                }
            }
            else
            {
                ISAAC_ASSERT_MSG(false, "Overlapping gaps are not allowed");
//                unsigned appliedGapLength = gap.getEndPos(true) - lastGapEndPos;
//                if (gap.isInsertion())
//                {
//                    appliedGapLength = std::min<unsigned>(appliedGapLength, basesLeft);
//                    basesLeft -= appliedGapLength;
//                }
//                realignedCigars_.push_back(Cigar::encode(appliedGapLength, gap.getOpCode()));
//                lastOperation = gap.getOpCode();
            }
        }
        ++currentGapIndex;
    }
    if(basesLeft)
    {
        ISAAC_ASSERT_MSG2(0 < basesLeft, "Spent more than readLength (%d left) bases on the CIGAR: %s",
                          basesLeft % alignment::Cigar::toString(&realignedCigars_.at(before), &realignedCigars_.back() + 1));
        ISAAC_ASSERT_MSG(realignedCigars_.capacity() > realignedCigars_.size(), "Realigned CIGAR buffer is out of capacity");
        realignedCigars_.push_back(Cigar::encode(basesLeft, Cigar::ALIGN));
    }

    if (fragment.rightClipped())
    {
        realignedCigars_.push_back(Cigar::encode(fragment.rightClipped(), Cigar::SOFT_CLIP));
    }

    index.pos_ = newBeginPos;
    index.cigarBegin_ = &realignedCigars_.at(before);
    index.cigarEnd_ = &realignedCigars_.back() + 1;
}

/**
 * \brief Perform full realignment discarding all the existing gaps
 */
void GapRealigner::realignFast(
    const reference::ReferencePosition binStartPos,
    const reference::ReferencePosition binEndPos,
    PackedFragmentBuffer::Index &index,
    PackedFragmentBuffer &dataBuffer)
{
    realign(binStartPos, binEndPos, index, dataBuffer);
}

/**
 * \brief Find the start position such that the base that would be the read base
 *        at pivotPos if read originally had no gaps, would still be at pivotPos
 *        given all gaps in this choice are applied
 */
 

bool GapRealigner::findStartPos(
    const unsigned short choice,
    const GapsRange &gaps,
    const reference::ReferencePosition binStartPos,
    const PackedFragmentBuffer::Index &index,
    const unsigned pivotGapIndex,
    const reference::ReferencePosition pivotPos,
    reference::ReferencePosition &ret)
{
    // Find the start position such that the base at pivotPos does not move when all
    // existing fragment gaps are removed
    using alignment::Cigar;
    reference::ReferencePosition lastGapEndPos = index.pos_;

    // initialize with the distance between the (possibly clipped) read start and pivotPos,
    // then offset by the gaps from the original CIGAR.
    long offset = std::max(0L, pivotPos - index.pos_);
//    ISAAC_THREAD_CERR << " restoring startPos offset=" << offset << " from pivot=" << pivotPos << "and original=" << index.pos_ << std::endl;
    BOOST_FOREACH(const uint32_t cigarOp, std::make_pair(index.cigarBegin_, index.cigarEnd_))
    {
        if (lastGapEndPos >= pivotPos)
        {
            break;
        }
        const Cigar::Component decoded = Cigar::decode(cigarOp);
        if (Cigar::ALIGN == decoded.second)
        {
            lastGapEndPos += decoded.first;
        }
        else if (Cigar::INSERT == decoded.second)
        {
            offset += decoded.first;
        }
        else if (Cigar::DELETE == decoded.second)
        {
            lastGapEndPos += decoded.first;
            offset -= decoded.first;
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
            offset += decoded.first;
        }
        else
        {
            ISAAC_ASSERT_MSG2(false, "Unexpected CIGAR operation %d in %0xd", decoded.second % cigarOp);
        }
    }

//    ISAAC_THREAD_CERR << " restored startPos offset=" << offset << " from pivot=" << pivotPos << "and original=" << index.pos_ << std::endl;

//    if (0 > offset)
//    {
        // pivot pos before the ungapped read start
//        return false;
//    }

    // Now shift found start position by the gaps included in the mask while keeping
    // the base at pivotPos steady
    unsigned gapIndex = pivotGapIndex - 1;
    reference::ReferencePosition overlapPos = pivotPos;
    unsigned basesLeft = offset;
    BOOST_REVERSE_FOREACH(const gapRealigner::Gap& gap, std::make_pair(gaps.first, gaps.first + pivotGapIndex))
    {
        if (choice & (1 << gapIndex))
        {
            if (gap.getEndPos(false) > overlapPos)
            {
//                ISAAC_THREAD_CERR << " overlapping gaps are not allowed " << overlapPos << "-" << gap << std::endl;
                // overlapping gaps are not allowed
                return false;
            }
            if (gap.isInsertion())
            {
                const unsigned insertionBases = std::min(basesLeft, gap.getLength());
                offset -= insertionBases;
                basesLeft -= insertionBases;
                if (!basesLeft)
                {
                    // insertions have eaten all the bases
                    break;
                }
            }
            else
            {
                //ISAAC_THREAD_CERR << " thinking of shift" << pivotPos << "-" << ret << std::endl;
                offset += gap.getLength();
                //ISAAC_THREAD_CERR << "shift" << pivotPos << "-" << ret << std::endl;

                overlapPos = gap.getBeginPos();
            }
        }
        --gapIndex;
    }

    if (binStartPos + offset > pivotPos)
    {
//        ISAAC_THREAD_CERR << " gap places read before bin start" << std::endl;
        // this combination of gaps will have the read start position moved before the binStartPos.
        // Don't realign this way.
        return false;
    }

    ret = pivotPos - offset;
//    ISAAC_THREAD_CERR << " new startPos offset=" << offset << " startPos=" << ret << " from pivot=" << pivotPos << "and original=" << index.pos_ << std::endl;

    return true;
}

unsigned short getTotalGapsLength(
    const uint32_t *cigarIterator,
    const uint32_t * const cigarEnd)
{
    using alignment::Cigar;
    unsigned short ret = 0;
    for (;cigarIterator != cigarEnd; ++cigarIterator)
    {
        Cigar::Component decoded = Cigar::decode(*cigarIterator);
        if (Cigar::ALIGN == decoded.second)
        {
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
        }
        else if (Cigar::INSERT == decoded.second)
        {
            ret += decoded.first;
        }
        else if (Cigar::DELETE == decoded.second)
        {
            ret += decoded.first;
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation");
        }
    }

    return ret;
}

/**
 * \brief Perform full realignment discarding all the existing gaps
 */
void GapRealigner::realign(
    const reference::ReferencePosition binStartPos,
    const reference::ReferencePosition binEndPos,
    PackedFragmentBuffer::Index &index,
    PackedFragmentBuffer &dataBuffer)
{
    using alignment::Cigar;
    io::FragmentAccessor &fragment = dataBuffer.getFragment(index);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign " << fragment);
    // don't bother to look at perfectly matching ones
    // TODO: we can't update template lengths if mate is in a different bin. If realignment requires
    // template update and mate is in a different bin, the realignment should not be done. At the
    // moment just not realigning bin-spanning pairs at all.
    if (!fragment.flags_.unmapped_ && fragment.editDistance_ && fragment.flags_.mateBinTheSame_ &&
        // CASAVA IndelFinder bam reader skips reads containing soft-clip operations. This becomes
        // lethal with singleton-shadow pairs as then it pairs the shadow with something else (following read I guess)
        // In cases when it pairs such a orphaned shadow with an end of a chimeric read it then produces a monster pair
        // that even the almighty ClusterMerger cannot swallow.
        // Avoid realigning singletons for now. TODO: don't realign singletons only if doing so produces soft-clipping
        !fragment.flags_.mateUnmapped_ &&
        // Normal alignments don't hit gaps by large. Zero-scored templates tend to pile up around single locations.
        // Although with gcc -O3 this passes unnoticed, it causes enormous time waste trying to realign heap of
        // misplaced reads against the gaps they happen to have. TODO: enable the line below.
        ((unsigned short)(-1) != fragment.alignmentScore_ || (unsigned short)(-1) != fragment.templateAlignmentScore_))

    {
        const std::vector<reference::Contig> &reference = contigList_.at(barcodeMetadataList_.at(fragment.barcode_).getReferenceIndex());

        // Mate realignment might have updated our fragment.fStrandPosition_. Make sure index.pos_ is up to date
        index.pos_ = fragment.fStrandPosition_;
        const unsigned sampleId = barcodeBamMapping_.getFileIndex(fragment.barcode_);
        RealignmentBounds bounds = extractRealignmentBounds(index);
        const GapsRange gaps = findGaps(sampleId, binStartPos, bounds.beginPos_, bounds.endPos_);
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Found Gaps " << gaps << " for bounds " << bounds);

        if (MAX_GAPS_AT_A_TIME < gaps.size())
        {
            ISAAC_THREAD_CERR_DEV_TRACE("GapRealigner::realign: Too many gaps (" << gaps.size() << "). Will not realign: " << fragment);
            // computing all combinations will simply take too long.
            return;
        }
        //ISAAC_ASSERT_MSG2(gaps.size() < MAX_GAPS_AT_A_TIME, "Too many gaps for realignment: %d", gaps.size());
        const unsigned short maxChoice = (1 << gaps.size()) - 1;
        unsigned bestEditDistance = fragment.editDistance_;
        unsigned bestMismatches = fragment.editDistance_ - getTotalGapsLength(index.cigarBegin_, index.cigarEnd_);
//        ISAAC_THREAD_CERR << "Initial bestEditDistance " << bestEditDistance <<
//            " bestMismatches " << bestMismatches << std::endl;

        reference::ReferencePosition bestStartPos = index.pos_;
        // 0 means none of the gaps apply. It also means the original alignment should be kept.
        unsigned short bestChoice = 0;
        for (unsigned short choice = 0; choice <= maxChoice; ++choice)
        {
            unsigned pivotGapIndex = 0;
            BOOST_FOREACH(const gapRealigner::Gap& pivotGap, std::make_pair(gaps.first, gaps.second))
            {
                if (choice & (1 << pivotGapIndex))
                {
                    //TODO: check binBorder overrun
                    reference::ReferencePosition newStarPos;

                    // verify case when anchoring occurs before pivot gap
                    if (pivotGap.getBeginPos() >= binStartPos)
                    {
//                        ISAAC_THREAD_CERR << "Testing choice " << int(choice) <<
//                            " pivot before " << pivotGap << std::endl;
                        if (findStartPos(choice, gaps, binStartPos, index, pivotGapIndex,
                                         pivotGap.getBeginPos(), newStarPos))
                        {
                            unsigned thisChoiceEditDistance = 0;
                            const unsigned thisChoiceMismatches =
                                verifyGapsChoice(choice, gaps, newStarPos, binStartPos,
                                                       index, fragment, bounds, reference,
                                                       thisChoiceEditDistance);
//                            ISAAC_THREAD_CERR << "Tested choice " << int(choice) << " ed=" <<
//                                thisChoiceEditDistance << " mm=" << thisChoiceMismatches << " pivot before " << pivotGap <<
//                                " new start pos " << newStarPos << std::endl;
                            if (thisChoiceMismatches < bestMismatches ||
                                (thisChoiceMismatches == bestMismatches &&
                                    thisChoiceEditDistance <= bestEditDistance))
                            {
//                                ISAAC_THREAD_CERR << "Elected choice " << int(choice) << " ed=" <<
//                                    thisChoiceEditDistance << " mm=" << thisChoiceMismatches << " pivot before " << pivotGap <<
//                                    " new start pos " << newStarPos << std::endl;
                                bestEditDistance = thisChoiceEditDistance;
                                bestChoice = choice;
                                bestStartPos = newStarPos;
                                bestMismatches = thisChoiceMismatches;
                            }
                        }
                    }

//                    ISAAC_THREAD_CERR << "Testing choice " << int(choice) <<
//                        " pivot after " << pivotGap << std::endl;

                    // verify case when anchoring occurs after pivot gap
                    if (findStartPos(choice, gaps, binStartPos, index, pivotGapIndex + 1,
                                     pivotGap.getEndPos(false), newStarPos))
                    {
                        unsigned thisChoiceEditDistance = 0;
                        const unsigned thisChoiceMismatches =
                            verifyGapsChoice(choice, gaps, newStarPos, binStartPos,
                                                   index, fragment, bounds, reference,
                                                   thisChoiceEditDistance);
//                        ISAAC_THREAD_CERR << "Tested choice " << int(choice) << " ed=" <<
//                            thisChoiceEditDistance << " mm=" << thisChoiceMismatches <<  " pivot after " << pivotGap <<
//                            " new start pos " << newStarPos << std::endl;
                        if (thisChoiceMismatches < bestMismatches ||
                            (thisChoiceMismatches == bestMismatches &&
                                thisChoiceEditDistance <= bestEditDistance))
                        {
//                            ISAAC_THREAD_CERR << "Elected choice " << int(choice) << " ed=" <<
//                                thisChoiceEditDistance << " mm=" << thisChoiceMismatches <<  " pivot after " << pivotGap <<
//                                " new start pos " << newStarPos << std::endl;

                            bestEditDistance = thisChoiceEditDistance;
                            bestChoice = choice;
                            bestStartPos = newStarPos;
                            bestMismatches = thisChoiceMismatches;
                        }
                    }
                }
                ++pivotGapIndex;
            }
        }

        if (bestChoice && binEndPos > bestStartPos)
        {
//            ISAAC_THREAD_CERR << "Applying choice for bin end pos:" << binEndPos << int(bestChoice) << " to gaps " << gaps << index << fragment << std::endl;
            PackedFragmentBuffer::Index tmp = index;
            tmp.pos_ = bestStartPos;
            applyChoice(bestChoice, gaps, binStartPos, tmp, fragment, dataBuffer);

//            ISAAC_THREAD_CERR << " before compactCigar=" << index << fragment << std::endl;
            const unsigned short oldEditDistance = fragment.editDistance_;
            fragment.editDistance_= bestEditDistance;
            if (compactCigar(binEndPos, tmp, fragment))
            {
                if (clipSemialigned_)
                {
                    // Note! this has to be called after compactCigar as otherwise the fragment.observedLength_ is incorrect
                    SemialignedEndsClipper clipper(realignedCigars_);
                    clipper.clip(reference, binEndPos, tmp, fragment);
                }

                index = tmp;
//                ISAAC_THREAD_CERR << " before updatePairDetails=" << index << fragment << std::endl;
                updatePairDetails(index, fragment, dataBuffer);
//                ISAAC_THREAD_CERR << "Applying choice done " << int(bestChoice) << " to gaps " << gaps << index << fragment << std::endl;
            }
            else
            {
                fragment.editDistance_ = oldEditDistance;
//                ISAAC_THREAD_CERR << "Ignoring choice" << int(bestChoice) << " to gaps " << gaps << index << fragment <<
//                    "due to cigar compacting having to move read into the next bin." << std::endl;
            }
        }
        else
        {
//            ISAAC_THREAD_CERR << "Ignoring choice" << int(bestChoice) << " to gaps " << gaps << index << fragment <<
//                "due to realigned read having to move to the next bin." << std::endl;
        }
    }
    else
    {
        ISAAC_THREAD_CERR_DEV_TRACE("GapRealigner::realign: Will not realign " << fragment);
    }
}


} // namespace build
} // namespace isaac
