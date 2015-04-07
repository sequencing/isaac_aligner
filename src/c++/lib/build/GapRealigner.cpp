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
#include "build/gapRealigner/OverlappingGapsFilter.hh"

#include "SemialignedEndsClipper.hh"

namespace isaac
{
namespace build
{

void RealignerGaps::reserve(const size_t gaps)
{
    deletionEndGroups_.reserve(gaps);
    gapGroups_.reserve(gaps);
}

void RealignerGaps::unreserve()
{
    gapRealigner::Gaps().swap(gapGroups_);
    gapRealigner::Gaps().swap(deletionEndGroups_);
}


inline bool orderByGapStartAndTypeLength(const gapRealigner::Gap &left, const gapRealigner::Gap &right)
{
    // ordering by signed length is required to ensure that intermixing of insertions and deletions does not prevent
    // duplicate gaps from being removed in finalizeGaps.
    return left.getBeginPos() < right.getBeginPos() ||
        (left.getBeginPos() == right.getBeginPos() && left.length_ < right.length_);
}

inline bool orderByDeletionGapEnd(const gapRealigner::Gap &left, const gapRealigner::Gap &right)
{
    return left.getDeletionEndPos() < right.getDeletionEndPos();
}

inline std::ostream & operator << (std::ostream &os, const GapRealigner::RealignmentBounds &fragmentGaps)
{
    return os << "FragmentGaps(" <<
        fragmentGaps.beginPos_ << "," <<
        fragmentGaps.firstGapStartPos_ << "," <<
        fragmentGaps.lastGapEndPos_ << "," <<
        fragmentGaps.endPos_ << ")";
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

void RealignerGaps::finalizeGaps()
{
    std::sort(gapGroups_.begin(), gapGroups_.end(), orderByGapStartAndTypeLength);
    gapGroups_.erase(std::unique(gapGroups_.begin(), gapGroups_.end()), gapGroups_.end());

    std::remove_copy_if(gapGroups_.begin(), gapGroups_.end(), std::back_inserter(deletionEndGroups_),
                        !boost::bind(&gapRealigner::Gap::isDeletion, _1));
    std::sort(deletionEndGroups_.begin(), deletionEndGroups_.end(), orderByDeletionGapEnd);
}

/**
 * \brief Find gaps given the position range.
 */
gapRealigner::GapsRange RealignerGaps::findGaps(
    const unsigned long clusterId,
    const reference::ReferencePosition binStartPos,
    const reference::ReferencePosition rangeBegin,
    const reference::ReferencePosition rangeEnd,
    gapRealigner::Gaps &foundGaps) const
{
    foundGaps.clear();

    //ISAAC_THREAD_CERR << "findGaps effective range: [" << earliestPossibleGapBegin << ";" << rangeEnd << ")" << std::endl;
    gapRealigner::GapsRange gapStarts;
//    ISAAC_THREAD_CERR_DEV_TRACE("findGaps all gaps: " << gapRealigner::GapsRange(sampleGaps.begin(), sampleGaps.end()));
    // the first one that begins on or after the rangeBegin
    gapStarts.first = std::lower_bound(gapGroups_.begin(), gapGroups_.end(),
                                 gapRealigner::Gap(rangeBegin, -1000000), orderByGapStartAndTypeLength);
    // the first one that ends on or after the rangeEnd
    gapStarts.second = std::lower_bound(gapStarts.first, gapGroups_.end(), gapRealigner::Gap(rangeEnd, 0), orderByGapStartAndTypeLength);

    gapRealigner::GapsRange gapEnds;
    gapEnds.first = std::lower_bound(deletionEndGroups_.begin(), deletionEndGroups_.end(),
                                     // gap length 1 is simply to prevent orderByDeletionGapEnd for failing an assertion
                                     gapRealigner::Gap(rangeBegin, 1), orderByDeletionGapEnd);
    gapEnds.second = std::lower_bound(gapEnds.first, deletionEndGroups_.end(),
                                      // gap length 1 is simply to prevent orderByDeletionGapEnd for failing an assertion
                                      gapRealigner::Gap(rangeEnd, 1), orderByDeletionGapEnd);

    if (foundGaps.capacity() < gapStarts.size() + gapEnds.size())
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, "GapRealigner::findGaps: Too many gaps (" << gapStarts.size() << "+" << gapEnds.size() << ")");
    }
    else
    {
        foundGaps.insert(foundGaps.end(), gapStarts.first, gapStarts.second);
        foundGaps.insert(foundGaps.end(), gapEnds.first, gapEnds.second);

        if (!gapEnds.empty() && !gapStarts.empty())
        {
            // if both ranges contain some gaps, we need to consolidate...
            std::sort(foundGaps.begin(), foundGaps.end(), orderByGapStartAndTypeLength);
            foundGaps.erase(std::unique(foundGaps.begin(), foundGaps.end()), foundGaps.end());
        }
    }

    return gapRealigner::GapsRange(foundGaps.begin(), foundGaps.end());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
                ret.lastGapEndPos_ = ret.firstGapStartPos_ = (ret.beginPos_ -= decoded.first);
            }
            else
            {
                ISAAC_ASSERT_MSG(index.cigarEnd_ == cigarIterator + 1,
                                 "At most two soft-clips are expected with second one being the last component of the cigar");
                ret.endPos_ += decoded.first;
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
            ret.endPos_ += decoded.first;
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
    const unsigned compareLength = std::min<unsigned>(length, std::distance(referenceBaseIt, contig.forward_.end()));
    unsigned mismatches = 0;
    BOOST_FOREACH(const unsigned char readBase, std::make_pair(basesIterator, basesIterator + compareLength))
    {
        mismatches += *referenceBaseIt != oligo::getUppercaseBaseFromBcl(readBase);
        ++referenceBaseIt;
    }
/*
    referenceBaseIt = contig.forward_.begin() + pos.getPosition();

    ISAAC_THREAD_CERR << mismatches << " mismatches " << compareLength << "compareLength " << pos <<
        " read '" << oligo::bclToString(basesIterator, compareLength) <<
        "' ref '" << std::string(referenceBaseIt, referenceBaseIt + compareLength) << "'" <<
        std::endl;
*/
    return mismatches;
}

/**
 * \brief Adjusts template length and mate fields
 *
 * \param beginPosShift change in the fragment alignment position. negative for insertions, positive for deletions
 * \param endPosShift   change in the fragment end alignment position. negative for insertions, positive for deletions
 */
void GapRealigner::updatePairDetails(
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
            ISAAC_ASSERT_MSG(mate.bamTlen_ == -oldBamTlen,
                             "Incoherent BAM template length between fragment and unmapped mate " << index <<
                             " " << fragment << " mate " << mate);
            mate.bamTlen_ = -fragment.bamTlen_;
            fragment.mateFStrandPosition_ = fragment.fStrandPosition_;
            mate.mateFStrandPosition_ = fragment.fStrandPosition_;
            mate.fStrandPosition_ = fragment.fStrandPosition_;
        }
        return;
    }

    io::FragmentAccessor &mate  = dataBuffer.getMate(index);

    ISAAC_ASSERT_MSG(mate.bamTlen_ == -fragment.bamTlen_,
                     "Incoherent BAM template length between fragment and mate " << index <<
                     " " << fragment <<
                     " mate " << mate);

    const reference::ReferencePosition fragmentBeginPos = fragment.fStrandPosition_;
    const reference::ReferencePosition fragmentEndPos = fragmentBeginPos + fragment.observedLength_;
    const reference::ReferencePosition mateBeginPos = fragment.mateFStrandPosition_;
    const reference::ReferencePosition mateEndPos = mateBeginPos + mate.observedLength_;

    fragment.bamTlen_ = io::FragmentAccessor::getTlen(fragmentBeginPos, fragmentEndPos, mateBeginPos, mateEndPos, fragment.flags_.firstRead_);
    mate.bamTlen_ = -fragment.bamTlen_;
    mate.mateFStrandPosition_ = fragment.fStrandPosition_;
    mate.flags_.properPair_ = fragment.flags_.properPair_ =
        alignment::TemplateLengthStatistics::Nominal ==
            barcodeTemplateLengthStatistics_.at(fragment.barcode_).checkModel(fragment, mate);
}

/**
 * \breif Collapses gaps on the ends of CIGAR into soft-clips
 *        This keeps CASAVA happy and probably in general is a right thing to do.
 *        Also adjusts fragment position and observed length.
 *
 * \return false, when the compacting the cigar would produce invalid fragment. Currently two cases are considered invalid:
 *          a) move the read past the binEndPos
 *          b) fragment gets completely soft-clipped
 *         If false is returned index and fragment are guaranteed to be unchanged.
 */
bool GapRealigner::compactCigar(
    const std::vector<reference::Contig> &reference,
    const reference::ReferencePosition binEndPos,
    PackedFragmentBuffer::Index &index,
    io::FragmentAccessor &fragment)
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, " Compacting " << fragment << index);

    ISAAC_ASSERT_MSG(alignment::Cigar::getReadLength(index.cigarBegin_, index.cigarEnd_) == fragment.readLength_,
                     "Broken CIGAR before compacting: " << alignment::Cigar::toString(index.cigarBegin_, index.cigarEnd_) <<
                     " " << index << " " << fragment);

    // at this point the fragment.editDistance_ is not adjusted for soft clipped bases
    using alignment::Cigar;
    const uint32_t *cigarIterator = index.cigarBegin_;
    unsigned softClipStart = 0;
    bool needCompacting = false;
    reference::ReferencePosition newPos = index.pos_;
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
        }
        else if (Cigar::DELETE == decoded.second)
        {
            needCompacting = true;
            if (binEndPos <= newPos + decoded.first)
            {
                return false;
            }
            newPos += decoded.first;
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation");
        }
    }

    if (index.cigarEnd_ == cigarIterator)
    {
        // fragment gets completely soft-clipped
        return false;
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
                             "At most two soft-clips are expected with second one being the last component of the cigar" << index);
            softClipEnd += decoded.first;
        }
        else if (Cigar::INSERT == decoded.second)
        {
            needCompacting = true;
            softClipEnd += decoded.first;
        }
        else if (Cigar::DELETE == decoded.second)
        {
            needCompacting = true;
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
        ISAAC_ASSERT_MSG(std::distance(cigarIterator, cigarBackwardsIterator + 1) < 200 &&
                         std::distance(cigarIterator, cigarBackwardsIterator + 1) > 0,
                         "Suspiciously long (" << std::distance(cigarIterator, cigarBackwardsIterator + 1) <<
                         ") CIGAR in the middle of compacting: " <<
                         " " << alignment::Cigar::toString(index.cigarBegin_, index.cigarEnd_) <<
                         " " << index << " " << fragment);

        ISAAC_ASSERT_MSG(alignment::Cigar::getReadLength(cigarIterator, cigarBackwardsIterator + 1) + softClipStart + softClipEnd == fragment.readLength_,
                         "Broken CIGAR in the middle of compacting: " << alignment::Cigar::toString(cigarIterator, cigarBackwardsIterator + 1) <<
                         " " << index << " " << fragment);
        ISAAC_ASSERT_MSG(realignedCigars_.capacity() >= (realignedCigars_.size() + std::distance(cigarIterator, cigarBackwardsIterator + 1)), "Realigned CIGAR buffer is out of capacity");
        realignedCigars_.insert(realignedCigars_.end(), cigarIterator, cigarBackwardsIterator + 1);
        if (softClipEnd)
        {
            ISAAC_ASSERT_MSG(realignedCigars_.capacity() > realignedCigars_.size(), "Realigned CIGAR buffer is out of capacity");
            realignedCigars_.push_back(Cigar::encode(softClipEnd, Cigar::SOFT_CLIP));
        }

        index.cigarBegin_ = &realignedCigars_.at(before);
        index.cigarEnd_ = &realignedCigars_.back() + 1;
        index.pos_ = newPos;
    }

    // recompute editDistance and observed length
    unsigned short newEditDistance = 0;
    const unsigned char *basesIterator = fragment.basesBegin() + softClipStart;
    std::vector<char>::const_iterator referenceIterator =
        reference.at(index.pos_.getContigId()).forward_.begin() + index.pos_.getPosition();

    reference::ReferencePosition newEndPos = index.pos_;
    for (;cigarIterator != cigarBackwardsIterator + 1; ++cigarIterator)
    {
        Cigar::Component decoded = Cigar::decode(*cigarIterator);
        if (Cigar::ALIGN == decoded.second)
        {
            newEditDistance += countMismatches(reference, basesIterator, newEndPos, decoded.first);
            newEndPos += decoded.first;
            basesIterator += decoded.first;
            referenceIterator += decoded.first;
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
            ISAAC_ASSERT_MSG(false, "At most two soft-clips are expected. Both at the ends of the CIGAR " << index << " " << fragment);
        }
        else if (Cigar::INSERT == decoded.second)
        {
            newEditDistance += decoded.first;
            basesIterator += decoded.first;
        }
        else if (Cigar::DELETE == decoded.second)
        {
            newEditDistance += decoded.first;
            newEndPos += decoded.first;
            referenceIterator += decoded.first;
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation");
        }
    }

    ISAAC_ASSERT_MSG(alignment::Cigar::getReadLength(index.cigarBegin_, index.cigarEnd_) == fragment.readLength_,
                      "Broken CIGAR after compacting: " << alignment::Cigar::toString(index.cigarBegin_, index.cigarEnd_) <<
                      " " << index << " " << fragment);

    fragment.editDistance_ = newEditDistance;
    fragment.fStrandPosition_ = index.pos_;
    fragment.observedLength_ = newEndPos - fragment.fStrandPosition_;

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, " Compacted " << fragment << index);

    return true;
}

/**
 * \brief bits in choice determine whether the corresponding gaps are on or off
 *
 * \return cost of the new choice or -1U if choice is inapplicable.
 */
GapRealigner::GapChoice GapRealigner::verifyGapsChoice(
    const unsigned short choice,
    const gapRealigner::GapsRange &gaps,
    const reference::ReferencePosition newBeginPos,
    const PackedFragmentBuffer::Index &index,
    const io::FragmentAccessor &fragment,
    const std::vector<reference::Contig> &reference)
{
    GapChoice ret;
    // keeping as int to allow debug checks for running into negative
    int basesLeft = fragment.readLength_;
    int leftClippedLeft = fragment.leftClipped();

//    newBeginPos += fragment.leftClipped();
//    basesLeft -= fragment.leftClipped();

    reference::ReferencePosition lastGapEndPos = newBeginPos;
    reference::ReferencePosition lastGapBeginPos; // initially set to an invalid position which would not match any gap pos
    unsigned currentGapIndex = 0;
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
                ret.cost_ = -1U;
                return ret;
            }
//            ISAAC_THREAD_CERR << " lastGapEndPos=" << lastGapEndPos << std::endl;

            if (gap.getBeginPos() < lastGapEndPos)
            {
                ret.cost_ = -1U;
                return ret;
                // Allowing overlapping deletions is tricky because it is hard to track back the
                // newBeginPos from the pivot see findStartPos.
            }

            if (gap.getBeginPos() == lastGapBeginPos)
            {
                // The only case where it makes sense to allow two or more gaps starting at the same
                // position is when we want to combine multiple insertions into a larger
                // one. Unfortunately, with enough gaps, it consumes the read into one single insertion...
                // Other cases:
                // deletion/deletion - is disallowed above
                // insertion/deletion - does not make sense (and cause trouble SAAC-253)
                ret.cost_ = -1U;
                return ret;
            }

//            ISAAC_THREAD_CERR << " lastGapEndPos=" << lastGapEndPos << " lastGapBeginPos=" << lastGapBeginPos <<
//                " gap.isInsertion()=" << gap.isInsertion() << std::endl;

            const int mappedBases = std::min<int>(basesLeft - fragment.rightClipped(), gap.getBeginPos() - lastGapEndPos);
//            ISAAC_THREAD_CERR << " mappedBases=" << mappedBases << " basesLeft=" << basesLeft << std::endl;

            const unsigned length = mappedBases - std::min(mappedBases, leftClippedLeft);
            const unsigned mm = countMismatches(reference,
                                                fragment.basesBegin() + (fragment.readLength_ - basesLeft) + leftClippedLeft,
                                                lastGapEndPos + leftClippedLeft, length);

//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "countMismatches: " << mm);
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "leftClippedLeft: " << leftClippedLeft);

            ret.mappedLength_ += length;
            ret.editDistance_ += mm;
            ret.mismatches_ += mm;
            ret.cost_ += mm * mismatchCost_;
            basesLeft -= mappedBases;
            leftClippedLeft -= std::min(leftClippedLeft, mappedBases);
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "leftClippedLeft: " << leftClippedLeft);
            unsigned clippedGapLength = 0;
            if (gap.isInsertion())
            {
                clippedGapLength = std::min<int>(basesLeft - fragment.rightClipped(), gap.getLength());
                // insertions reduce read length
                basesLeft -= clippedGapLength;
                leftClippedLeft -= std::min<int>(leftClippedLeft, gap.getLength());
            }
            else
            {
                clippedGapLength = leftClippedLeft ? 0 : gap.getLength();
            }
            ret.editDistance_ += clippedGapLength;
            ret.cost_ += clippedGapLength ? (gapOpenCost_ + (clippedGapLength - 1) * gapExtendCost_) : 0;
            lastGapEndPos = gap.getEndPos(false);
            lastGapBeginPos = gap.getBeginPos();

            if (basesLeft == leftClippedLeft + fragment.rightClipped())
            {
                break;
            }
            ISAAC_ASSERT_MSG(basesLeft > leftClippedLeft + fragment.rightClipped(), "Was not supposed to run into the clipping");

        }
        ++currentGapIndex;
    }

    if(basesLeft > leftClippedLeft + fragment.rightClipped())
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "leftClippedLeft: " << leftClippedLeft);

        const unsigned length = basesLeft - std::min<unsigned>(basesLeft, leftClippedLeft) - fragment.rightClipped();

        const reference::ReferencePosition firstUnclippedPos = lastGapEndPos + leftClippedLeft;
        if (firstUnclippedPos.getPosition() > reference.at(firstUnclippedPos.getContigId()).forward_.size())
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, ret << " gap pushes part of the read outside the reference " << firstUnclippedPos << " " << basesLeft);
            ret.cost_ = -1U;
            return ret;
        }
        else
        {
            const unsigned mm = countMismatches(
                reference, fragment.basesBegin() + (fragment.readLength_ - basesLeft) + leftClippedLeft, firstUnclippedPos, length);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "final countMismatches: " << mm);
            ret.mappedLength_ += length;
            ret.editDistance_ += mm;
            ret.mismatches_ += mm;
            ret.cost_ += mm * mismatchCost_;
        }
    }
    else
    {
        ISAAC_ASSERT_MSG(leftClippedLeft + fragment.rightClipped() == basesLeft, "Spent more than readLength. basesLeft: " << basesLeft <<
            " choice: " << int(choice) <<
            " " << index <<
            " " << fragment <<
            ", newBeginPos " << newBeginPos <<
            " lastGapEndPos " << lastGapEndPos <<
            " leftClippedLeft " << leftClippedLeft <<
            " gaps: " << gaps);
    }

    return ret;

}

/**
 * \return false if choice cannot be applied. Currently this can happen due to left-side clipping having to move the
 *         read into the next bin
 *
 * \postcondition When false is returned, all state variables and writable inputs retain their original values.
 */
bool GapRealigner::applyChoice(
    const unsigned short choice,
    const gapRealigner::GapsRange &gaps,
    const reference::ReferencePosition binEndPos,
    const reference::ReferencePosition contigEndPos,
    PackedFragmentBuffer::Index &index,
    const io::FragmentAccessor &fragment)
{
//    ISAAC_THREAD_CERR << "GapRealigner::applyChoice index=" << index << std::endl;
    reference::ReferencePosition newBeginPos = index.pos_;
    using alignment::Cigar;

    const size_t before = realignedCigars_.size();

    int basesLeft = fragment.readLength_;

    int leftClippedLeft = fragment.leftClipped();
    // since insertions overlapped by left clipping don't move the alignment position, we need to count the overlaps for the final alignment position adjustment
    int leftClippedInsertionBases = 0;

    if (fragment.leftClipped())
    {
        Cigar::Component decoded = Cigar::decode(*fragment.cigarBegin());
        ISAAC_ASSERT_MSG(Cigar::SOFT_CLIP == decoded.second, "Original CIGAR is expected to have soft clip at the start " << fragment);
        ISAAC_ASSERT_MSG(fragment.leftClipped() <= decoded.first, "Original CIGAR soft clip at the start is shorter than left-clipped bases");
        realignedCigars_.push_back(Cigar::encode(fragment.leftClipped(), Cigar::SOFT_CLIP));
//        newBeginPos += fragment.leftClipped();
//        basesLeft -= fragment.leftClipped();
    }

    if (fragment.rightClipped())
    {
        Cigar::Component decoded = Cigar::decode(*(fragment.cigarEnd() - 1));
        ISAAC_ASSERT_MSG(Cigar::SOFT_CLIP == decoded.second, "Original CIGAR is expected to have soft clip at the end " << fragment);
        ISAAC_ASSERT_MSG(fragment.rightClipped() <= decoded.first, "Original CIGAR soft clip at the end is shorter than right-clipped bases");
//        basesLeft -= fragment.rightClipped();
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
                const int mappedBases = std::min<int>(basesLeft - fragment.rightClipped(), gapClippedBeginPos - lastGapEndPos);
                const unsigned softClippedMappedLength = mappedBases - std::min(mappedBases, leftClippedLeft);
//                ISAAC_THREAD_CERR << "lastGapEndPos=" << lastGapEndPos << std::endl;
//                ISAAC_THREAD_CERR << "mappedBases=" << mappedBases << std::endl;
//                ISAAC_THREAD_CERR << "softClippedMappedLength=" << softClippedMappedLength << std::endl;

                if (softClippedMappedLength)
                {
                    // avoid 0M in CIGAR
                    ISAAC_ASSERT_MSG(realignedCigars_.capacity() > realignedCigars_.size(), "Realigned CIGAR buffer is out of capacity");
                    realignedCigars_.push_back(Cigar::encode(softClippedMappedLength, Cigar::ALIGN));
                }
                basesLeft -= mappedBases;
                leftClippedLeft -= std::min(mappedBases, leftClippedLeft);
//                ISAAC_THREAD_CERR << "mappedBases=" << mappedBases << std::endl;

                if (gap.isInsertion())
                {
                    const int clippedGapLength =
                        std::min<int>(basesLeft - fragment.rightClipped(), gap.getEndPos(true) - gapClippedBeginPos);
                    const int softClippedGapLength = clippedGapLength - std::min(clippedGapLength, leftClippedLeft);
//                    ISAAC_THREAD_CERR << "clippedGapLength=" << clippedGapLength << std::endl;
                    if (softClippedGapLength)
                    {
                        if (alignment::Cigar::INSERT == lastOperation && !mappedBases)
                        {
                            //GATK does not like 2I1I type cigars
                            realignedCigars_.back() =
                                Cigar::encode(Cigar::decode(realignedCigars_.back()).first + softClippedGapLength,
                                              alignment::Cigar::INSERT);
                        }
                        else
                        {
                            realignedCigars_.push_back(Cigar::encode(softClippedGapLength, gap.getOpCode()));
                            ISAAC_ASSERT_MSG(realignedCigars_.capacity() > realignedCigars_.size(), "Realigned CIGAR buffer is out of capacity");
                            lastOperation = alignment::Cigar::INSERT;
                        }
                    }

                    basesLeft -= clippedGapLength;
                    lastGapEndPos = gapClippedBeginPos;
                    leftClippedLeft -= std::min(clippedGapLength, leftClippedLeft);
                    leftClippedInsertionBases += clippedGapLength - softClippedGapLength;

//                    ISAAC_THREAD_CERR << "2nd lastGapEndPos=" << lastGapEndPos << std::endl;
                }
                else
                {
                    // if left-side alignment-indepndent soft-clipping is in place, the deletions that we put in will simply move the alignment position forward in compactCigar
                    const int clippedGapLength = gap.getEndPos(true) - gapClippedBeginPos;
//                    ISAAC_THREAD_CERR << "clippedGapLength=" << clippedGapLength << std::endl;
                    if (!leftClippedLeft)
                    {
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
                    }
                    else
                    {
                        newBeginPos += clippedGapLength;
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

            if (basesLeft == leftClippedLeft + fragment.rightClipped())
            {
                break;
            }
            ISAAC_ASSERT_MSG(basesLeft > leftClippedLeft + fragment.rightClipped(), "Was not supposed to run into the clipping");
        }
        ++currentGapIndex;
    }
    if(basesLeft > leftClippedLeft + fragment.rightClipped())
    {
        ISAAC_ASSERT_MSG(realignedCigars_.capacity() > realignedCigars_.size(), "Realigned CIGAR buffer is out of capacity");
        const int basesToTheEndOfContig = contigEndPos - lastGapEndPos - leftClippedLeft;
        const int mappedBases = std::min(basesToTheEndOfContig, basesLeft - leftClippedLeft - fragment.rightClipped());
        if (mappedBases)
        {
            // avoid 0M in CIGAR
            realignedCigars_.push_back(Cigar::encode(mappedBases, Cigar::ALIGN));
        }
        basesLeft -= leftClippedLeft + mappedBases;
        leftClippedLeft = 0;
    }
    else
    {
        ISAAC_ASSERT_MSG(fragment.rightClipped() + leftClippedLeft == basesLeft, "Spent more than readLength (" << basesLeft <<
            " left) bases on the CIGAR: " << alignment::Cigar::toString(&realignedCigars_.at(before), &realignedCigars_.back() + 1));
    }


    if (basesLeft)
    {
        realignedCigars_.push_back(Cigar::encode(basesLeft, Cigar::SOFT_CLIP));
    }

    newBeginPos += fragment.leftClipped() - leftClippedInsertionBases;
    if (newBeginPos >= binEndPos)
    {
        // all things considered, we can't apply this choice as it would place read into the next bin.
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "applying would move start " << newBeginPos << " past binEndPos " << binEndPos << index);
        realignedCigars_.resize(before);
        return false;
    }
    index.pos_ = newBeginPos;
    index.cigarBegin_ = &realignedCigars_.at(before);
    index.cigarEnd_ = &realignedCigars_.back() + 1;
    return true;
}

/**
 * \brief Find the start position such that the base that would be the read base
 *        at pivotPos if read originally had no gaps, would still be at pivotPos
 *        given all gaps in this choice are applied
 */
 

bool GapRealigner::findStartPos(
    const unsigned short choice,
    const gapRealigner::GapsRange &gaps,
    const reference::ReferencePosition binStartPos,
    const reference::ReferencePosition binEndPos,
    const PackedFragmentBuffer::Index &index,
    const unsigned pivotGapIndex,
    const reference::ReferencePosition pivotPos,
    reference::ReferencePosition &ret)
{
    // Find the start position such that the base at pivotPos does not move when all
    // existing fragment gaps are removed
    using alignment::Cigar;
    reference::ReferencePosition lastGapEndPos = index.getUnclippedPosition();

    // initialize with the distance between the (possibly clipped) read start and pivotPos,
    // then offset by the gaps and soft clipping from the original CIGAR.
    long offset = pivotPos - index.pos_;
//    ISAAC_THREAD_CERR << " restoring startPos offset=" << offset << " from pivot=" << pivotPos << "and original=" << index.pos_ << std::endl;
    for(const uint32_t *it = index.cigarBegin_; index.cigarEnd_ != it; ++it)
    {
        const uint32_t cigarOp = *it;
        if (lastGapEndPos > pivotPos)
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
            if (lastGapEndPos > pivotPos)
            {
                // existing deletion overlaps pivot position like this:
                //"A-----C",
                //"AGATCAG",
                //"   **");
                return false;
            }
            offset -= decoded.first;
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
            if(index.cigarBegin_ == it)
            {
                // soft clip at the start eats bases just like an insertion, but also moves the position same way as deletion does
                offset += decoded.first;
            }
            lastGapEndPos += decoded.first;
            // soft clip at the end is treated same way as mapped bases
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation " << decoded.second  << " in " << cigarOp);
        }
    }

//    ISAAC_THREAD_CERR << " restored startPos offset=" << offset << " from pivot=" << pivotPos << "and original=" << index.pos_ << std::endl;

    if (0 > offset)
    {
//        ISAAC_THREAD_CERR << " pivot pos before the ungapped read start " << offset << std::endl;
        // pivot pos before the ungapped read start
        return false;
    }

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
                ISAAC_THREAD_CERR << " overlapping gaps are not allowed " << overlapPos << "-" << gap << std::endl;
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

    if (pivotPos - offset >= binEndPos)
    {
//         ISAAC_THREAD_CERR << " gap places read after bin end" << std::endl;
         // this combination of gaps will have the read start position moved at or after the binEndPos.
         // Don't realign this way.
         return false;
    }

    ret = pivotPos - offset;
//    ISAAC_THREAD_CERR << " new startPos offset=" << offset << " startPos=" << ret << " from pivot=" << pivotPos << "and original=" << index.pos_ << std::endl;

    return true;
}

static unsigned short getTotalGapsLength(
    const uint32_t *cigarIterator,
    const uint32_t * const cigarEnd,
    unsigned &gapsCount,
    unsigned &mappedLength)
{
    gapsCount = 0;
    using alignment::Cigar;
    unsigned short ret = 0;
    for (;cigarIterator != cigarEnd; ++cigarIterator)
    {
        Cigar::Component decoded = Cigar::decode(*cigarIterator);
        if (Cigar::ALIGN == decoded.second)
        {
            mappedLength += decoded.first;
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
        }
        else if (Cigar::INSERT == decoded.second)
        {
            ret += decoded.first;
            ++gapsCount;
        }
        else if (Cigar::DELETE == decoded.second)
        {
            ret += decoded.first;
            ++gapsCount;
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation");
        }
    }

    return ret;
}

inline int calculateMismatchesPercent(
    unsigned mismatches,
    unsigned mappedLength)
{
//    ISAAC_THREAD_CERR << "calculateMismatchesPercent mismatches=" << mismatches << " mappedLength=" << mappedLength << std::endl;
    return mismatches * 100 / mappedLength;
}

unsigned GapRealigner::getAlignmentCost(
    const io::FragmentAccessor &fragment,
    const PackedFragmentBuffer::Index &index,
    unsigned &editDistance,
    int &mismatchesPercent) const
{
    unsigned gapsCount = 0;
    unsigned mappedLength = 0;
    const unsigned totalGapsLength = getTotalGapsLength(index.cigarBegin_, index.cigarEnd_, gapsCount, mappedLength);
    editDistance = fragment.editDistance_;
    const unsigned mismatches = fragment.editDistance_ - totalGapsLength;
//    ISAAC_THREAD_CERR << "getAlignmentCost " << fragment << std::endl;
    mismatchesPercent = calculateMismatchesPercent(mismatches, mappedLength);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Initial mm:" << mismatches << " gaps:" << gapsCount << " gapslength:" << totalGapsLength);
    return mismatches * mismatchCost_ + gapsCount * gapOpenCost_ + gapExtendCost_ * (totalGapsLength - gapsCount);
}

bool GapRealigner::isBetterChoice(
    const GapChoice &choice,
    const int originalMismatchesPercent,
    const unsigned bestCost,
    const unsigned bestEditDistance) const
{
//    ISAAC_THREAD_CERR << "isBetterChoice " << choice << std::endl;
    return
        choice.mappedLength_ &&
        (choice.cost_ < bestCost || (choice.cost_ == bestCost && choice.editDistance_ < bestEditDistance)) &&
        calculateMismatchesPercent(choice.mismatches_, choice.mappedLength_) <= originalMismatchesPercent;
}

/**
 * \brief Perform full realignment discarding all the existing gaps
 */
void GapRealigner::realign(
    RealignerGaps &realignerGaps,
    const reference::ReferencePosition binStartPos,
    reference::ReferencePosition binEndPos,
    PackedFragmentBuffer::Index &index,
    io::FragmentAccessor &fragment,
    PackedFragmentBuffer &dataBuffer)
{
    std::size_t bufferSizeBeforeRealignment = realignedCigars_.size();
    bool makesSenseToTryAgain = false;
//    bool firstAttempt = true;
//    lastAttemptGaps_.clear();

    if (fragment.flags_.unmapped_)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Will not realign unmapped " << fragment);
        return;
    }

    do
    {
        makesSenseToTryAgain = false;
        using alignment::Cigar;
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign " << index << fragment);
        const std::vector<reference::Contig> &reference = contigList_.at(barcodeMetadataList_.at(fragment.barcode_).getReferenceIndex());
        binEndPos = reference::ReferencePosition(
            binEndPos.getContigId(), std::min<unsigned long>(binEndPos.getPosition(), reference.at(binEndPos.getContigId()).getLength()));
        // don't bother to look at perfectly matching ones
        // TODO: we can't update template lengths if mate is in a different bin. If realignment requires
        // template update and mate is in a different bin, the realignment should not be done. At the
        // moment just not realigning bin-spanning pairs at all.
        if (fragment.editDistance_ &&
            // single-ended are ok to realign all the time
            (!fragment.flags_.paired_ ||
                // CASAVA IndelFinder bam reader skips reads containing soft-clip operations. This becomes
                // lethal with singleton-shadow pairs as then it pairs the shadow with something else (following read I guess)
                // In cases when it pairs such a orphaned shadow with an end of a chimeric read it then produces a monster pair
                // that even the almighty ClusterMerger cannot swallow.
                // Avoid realigning singletons for now. TODO: don't realign singletons only if doing so produces soft-clipping
                (!fragment.flags_.mateUnmapped_ &&
                    binStartPos <= fragment.mateFStrandPosition_ && binEndPos > fragment.mateFStrandPosition_)) &&
            // Normal alignments don't hit gaps by large. Zero-scored templates tend to pile up around single locations.
            // Although with gcc -O3 this passes unnoticed, it causes enormous time waste trying to realign heap of
            // misplaced reads against the gaps they happen to have.
            (realignDodgyFragments_ || fragment.DODGY_ALIGNMENT_SCORE != fragment.alignmentScore_ || fragment.DODGY_ALIGNMENT_SCORE != fragment.templateAlignmentScore_) &&
            // reference-clipped reads are difficult to compute because ReferencePosition is not allowed to go negative.
            (index.pos_.getPosition() >= index.getBeginClippedLength()))

        {

            // Mate realignment might have updated our fragment.fStrandPosition_. Make sure index.pos_ is up to date
            index.pos_ = fragment.fStrandPosition_;
            RealignmentBounds bounds = extractRealignmentBounds(index);
            const gapRealigner::GapsRange gaps = realignerGaps.findGaps(fragment.clusterId_, binStartPos, bounds.beginPos_, bounds.endPos_, currentAttemptGaps_);
//            if (!firstAttempt && lastAttemptGaps_.size() == currentAttemptGaps_.size() &&
//                std::equal(lastAttemptGaps_.begin(), lastAttemptGaps_.end(), currentAttemptGaps_.begin()))
//            {
//                // no new gaps found. stop trying.
//                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "no new gaps found. stop trying ");
//                break;
//            }
//            firstAttempt = false;
//            lastAttemptGaps_ = currentAttemptGaps_;

            if (!realignGapsVigorously_ && MAX_GAPS_AT_A_TIME < gaps.size())
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Too many gaps (" << gaps.size() << "). " << fragment);
                break;
            }

            const gapRealigner::OverlappingGapsFilter overlappingGapsFilter(gaps);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Found Gaps " << gaps << " for bounds " << bounds);

            unsigned bestEditDistance = 0;
            int originalMismatchesPercent = 0;
            unsigned bestCost = getAlignmentCost(fragment, index, bestEditDistance, originalMismatchesPercent);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Initial bestEditDistance " << bestEditDistance << " bestCost " << bestCost);

            reference::ReferencePosition bestStartPos = index.pos_;
            // 0 means none of the gaps apply. It also means the original alignment should be kept.
            unsigned bestChoice = 0;
            unsigned evaluatedSoFar = 0;
            for (unsigned choice = 0; (choice = overlappingGapsFilter.next(choice));)
            {
                if (((1 << MAX_GAPS_AT_A_TIME) - 1) < evaluatedSoFar++)
                {
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Too many gaps (" << evaluatedSoFar << " checked so far). " << fragment);
                    // We've spent too much time already. Just go with what we've got.
                    break;
                }

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
                            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Testing choice " << int(choice) << " pivot before " << pivotGap);
                            if (findStartPos(choice, gaps, binStartPos, binEndPos, index, pivotGapIndex,
                                             pivotGap.getBeginPos(), newStarPos))
                            {
                                const GapChoice thisChoice =
                                    verifyGapsChoice(choice, gaps, newStarPos, index, fragment, reference);
//                                ISAAC_THREAD_CERR << "Tested choice " << int(choice) << " ed=" <<
//                                    thisChoice.editDistance_ << " cost=" << thisChoice.cost_ << " pivot before " << pivotGap <<
//                                    " new start pos " << newStarPos << std::endl;
                                if (isBetterChoice(thisChoice, originalMismatchesPercent, bestCost, bestEditDistance))
                                {
                                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Elected choice " << int(choice) << " ed=" <<
                                        thisChoice.editDistance_ << " cost/best=" << thisChoice.cost_ << "/" << bestCost << " pivot before " << pivotGap <<
                                        " new start pos " << newStarPos);
                                    bestEditDistance = thisChoice.editDistance_;
                                    bestChoice = choice;
                                    bestStartPos = newStarPos;
                                    bestCost = thisChoice.cost_;
                                }
                                else
                                {
                                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Failed choice " << int(choice) << " ed=" <<
                                        thisChoice.editDistance_ << " cost/best=" << thisChoice.cost_ << "/" << bestCost << " pivot before " << pivotGap <<
                                        " new start pos " << newStarPos);
                                }
                            }
                        }

                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Testing choice " << int(choice) << " pivot after " << pivotGap);

                        // verify case when anchoring occurs after pivot gap
                        if (findStartPos(choice, gaps, binStartPos, binEndPos, index, pivotGapIndex + 1,
                                         pivotGap.getEndPos(false), newStarPos))
                        {
                            const GapChoice thisChoice =
                                verifyGapsChoice(choice, gaps, newStarPos, index, fragment, reference);
    //                        ISAAC_THREAD_CERR << "Tested choice " << int(choice) << " ed=" <<
    //                            thisChoiceEditDistance << " cost=" << thisChoiceCost <<  " pivot after " << pivotGap <<
    //                            " new start pos " << newStarPos << std::endl;
                            if (isBetterChoice(thisChoice, originalMismatchesPercent, bestCost, bestEditDistance))
                            {
                                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Elected choice " << int(choice) << " ed=" <<
                                    thisChoice.editDistance_ << " cost/best=" << thisChoice.cost_ << "/" << bestCost <<  " pivot after " << pivotGap <<
                                    " new start pos " << newStarPos);

                                bestEditDistance = thisChoice.editDistance_;
                                bestChoice = choice;
                                bestStartPos = newStarPos;
                                bestCost = thisChoice.cost_;
                            }
                            else
                            {
                                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Failed choice " << int(choice) << " ed=" <<
                                    thisChoice.editDistance_ << " cost/best=" << thisChoice.cost_ << "/" << bestCost <<  " pivot after " << pivotGap <<
                                    " new start pos " << newStarPos);
                            }
                        }
                    }
                    ++pivotGapIndex;
                }
            }

            if (bestChoice && binEndPos > bestStartPos)
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Applying choice for bin end pos:" << binEndPos <<
                                                       " " << int(bestChoice) << " to gaps " << gaps << index << fragment);
                PackedFragmentBuffer::Index tmp = index;
                tmp.pos_ = bestStartPos;

                const reference::ReferencePosition contigEndPos(binEndPos.getContigId(), reference.at(binEndPos.getContigId()).forward_.size());
                if (applyChoice(bestChoice, gaps, binEndPos, contigEndPos, tmp, fragment))
                {
        //            ISAAC_THREAD_CERR << " before compactCigar=" << index << fragment << std::endl;
                    if (compactCigar(reference, binEndPos, tmp, fragment))
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
                        makesSenseToTryAgain = realignGapsVigorously_;
                    }
                    else
                    {
        //                ISAAC_THREAD_CERR << "Ignoring choice" << int(bestChoice) << " to gaps " << gaps << index << fragment <<
        //                    "due to cigar compacting having to move read into the next bin." << std::endl;
                    }
                }
                else
                {
//                    ISAAC_THREAD_CERR << "Ignoring choice" << int(bestChoice) << " to gaps " << gaps << index << fragment <<
//                        "due to applyChoice having to move read into the next bin." << std::endl;
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
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Will not realign " << fragment);
        }
    } while(makesSenseToTryAgain);

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "before compactRealignedCigarBuffer:" << index);
    compactRealignedCigarBuffer(bufferSizeBeforeRealignment, index);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "after compactRealignedCigarBuffer:" << index);
}

void GapRealigner::compactRealignedCigarBuffer(
    const std::size_t bufferSizeBeforeRealignment,
    PackedFragmentBuffer::Index &index)
{
    if (realignedCigars_.size() != bufferSizeBeforeRealignment)
    {
        const std::size_t cigarLength = std::distance(index.cigarBegin_, index.cigarEnd_);
        const std::size_t expectedBufferSize = bufferSizeBeforeRealignment + cigarLength;
        if (expectedBufferSize  != realignedCigars_.size())
        {
            std::copy(index.cigarBegin_, index.cigarEnd_, realignedCigars_.begin() + bufferSizeBeforeRealignment);
            index.cigarBegin_ = &*(realignedCigars_.begin() + bufferSizeBeforeRealignment);
            index.cigarEnd_ = index.cigarBegin_ + cigarLength;
            realignedCigars_.resize(expectedBufferSize);
        }
    }
}


} // namespace build
} // namespace isaac
