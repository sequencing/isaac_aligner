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
 ** \file FragmentSequencingAdapterClipper.cpp
 **
 ** \brief See FragmentSequencingAdapterClipper.hh
 ** 
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "alignment/Alignment.hh"
#include "alignment/FragmentMetadata.hh"
#include "alignment/matchSelector/FragmentSequencingAdapterClipper.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

/**
 * \brief Adjusts the sequence iterators to stay within the reference.
 *
 * \return new fragment.position to point at the first not clipped base.
 */
long clipReference(
    const FragmentMetadata &fragmentMetadata,
    std::vector<char>::const_iterator &sequenceBegin,
    std::vector<char>::const_iterator &sequenceEnd,
    const long referenceSize)
{
    const long referenceLeft = referenceSize - fragmentMetadata.position;
    if (referenceLeft < std::distance(sequenceBegin, sequenceEnd))
    {
        sequenceEnd = sequenceBegin + referenceLeft;
    }

    if (0 > fragmentMetadata.position)
    {
        sequenceBegin -= fragmentMetadata.position;
        return 0L;
    }
    return fragmentMetadata.position;
}

unsigned countMatches(
    std::vector<char>::const_iterator sequenceBegin,
    const std::vector<char>::const_iterator sequenceEnd,
    std::vector<char>::const_iterator referenceBegin,
    const std::vector<char>::const_iterator referenceEnd)
{
    unsigned ret = 0;
    for (;sequenceEnd != sequenceBegin && referenceEnd != referenceBegin;
        ++sequenceBegin, ++ referenceBegin)
    {
        ret += isMatch(*sequenceBegin, *referenceBegin);
    }
    return ret;
}

/**
 * \return Number of !isMatch bases on the sequence/reference overlap
 */
unsigned countMismatches(
    std::vector<char>::const_iterator sequenceBegin,
    const std::vector<char>::const_iterator sequenceEnd,
    std::vector<char>::const_iterator referenceBegin,
    const std::vector<char>::const_iterator referenceEnd)
{
    return std::min(std::distance(sequenceBegin, sequenceEnd),
                    std::distance(referenceBegin, referenceEnd)) -
        countMatches(sequenceBegin, sequenceEnd, referenceBegin, referenceEnd);
}

/**
 * \return integer value representing % of !isMatch bases on the sequence/reference overlap
 */
unsigned percentMismatches(
    std::vector<char>::const_iterator sequenceBegin,
    const std::vector<char>::const_iterator sequenceEnd,
    std::vector<char>::const_iterator referenceBegin,
    const std::vector<char>::const_iterator referenceEnd)
{
    const unsigned overlapLength = std::min(std::distance(sequenceBegin, sequenceEnd),
                                            std::distance(referenceBegin, referenceEnd));
    return countMismatches(sequenceBegin, sequenceEnd, referenceBegin, referenceEnd) * 100 / overlapLength;
}


const std::pair<std::vector<char>::const_iterator, std::vector<char>::const_iterator>
findSequencingAdapter(
    std::vector<char>::const_iterator sequenceBegin,
    std::vector<char>::const_iterator sequenceEnd,
    const std::vector<char>::const_iterator referenceBegin,
    const std::vector<char>::const_iterator referenceEnd,
    const matchSelector::SequencingAdapter &adapter)
{
    for (std::vector<char>::const_iterator currentBase(sequenceBegin), currentReference(referenceBegin);
        sequenceEnd != currentBase; ++currentReference, ++currentBase)
    {
        if (!isMatch(*currentBase, *currentReference))
        {
            const std::pair<std::vector<char>::const_iterator, std::vector<char>::const_iterator> adapterMatchRange =
                adapter.getMatchRange(sequenceBegin, sequenceEnd, currentBase);
            if (adapterMatchRange.first != adapterMatchRange.second)
            {
                return adapterMatchRange;
            }
        }
    }
    return std::make_pair(sequenceBegin, sequenceBegin);
}

/**
 * \brief if the adapter sequence has not been checked for the alignment strand,
 *        searches for the adapter and prevents further searches for adapter on this strand
 */
void FragmentSequencingAdapterClipper::checkInitStrand(
    const FragmentMetadata &fragmentMetadata,
    const reference::Contig &contig)
{
    const bool reverse = fragmentMetadata.reverse;

    if (!strandAdapters_.strandRange_[reverse].initialized_)
    {
        std::vector<char>::const_iterator &adapterRangeBegin = strandAdapters_.strandRange_[reverse].adapterRangeBegin_;
        std::vector<char>::const_iterator &adapterRangeEnd = strandAdapters_.strandRange_[reverse].adapterRangeEnd_;

        const std::vector<char> &reference = contig.forward_;

        const Read &read = fragmentMetadata.getRead();

        const std::vector<char> &sequence = reverse ? read.getReverseSequence() : read.getForwardSequence();

        std::vector<char>::const_iterator sequenceBegin = sequence.begin();
        std::vector<char>::const_iterator sequenceEnd = sequence.end();

        long newFragmentPos =
            clipReference(fragmentMetadata, sequenceBegin, sequenceEnd, reference.size());

        adapterRangeBegin = sequenceEnd;
        adapterRangeEnd = sequenceBegin;

        BOOST_FOREACH(const matchSelector::SequencingAdapter &adapter, sequencingAdapters_)
        {
            if (adapter.isStrandCompatible(fragmentMetadata.isReverse()))
            {
                const std::pair<std::vector<char>::const_iterator, std::vector<char>::const_iterator> adapterMatchRange =
                    findSequencingAdapter(adapterRangeEnd, sequenceEnd,
                                          reference.begin() + newFragmentPos + std::distance(sequenceBegin, adapterRangeEnd),
                                          reference.end(), adapter);
                if (adapterMatchRange.first != adapterMatchRange.second)
                {
                    adapterRangeBegin = std::min(adapterMatchRange.first, adapterRangeBegin);
                    adapterRangeEnd = std::max(adapterMatchRange.second, adapterRangeEnd);
                    ISAAC_THREAD_CERR_DEV_TRACE("FragmentSequencingAdapterClipper::checkInitStrand found: " <<
                                                std::string(adapterRangeBegin, adapterRangeEnd) << " reverse: " << reverse);
                }
            }
        }
        strandAdapters_.strandRange_[reverse].initialized_ = true;
        strandAdapters_.strandRange_[reverse].empty_ = sequenceBegin == adapterRangeEnd;
    }
}

bool FragmentSequencingAdapterClipper::decideWhichSideToClip(
    const reference::Contig &contig,
    const long contigPosition,
    const std::vector<char>::const_iterator sequenceBegin,
    const std::vector<char>::const_iterator sequenceEnd,
    const SequencingAdapterRange &adapterRange,
    bool &clipBackwards)
{
    // try preserve the longest useful sequence
    const unsigned backwardsClipped = std::distance(sequenceBegin, adapterRange.adapterRangeBegin_);
    const unsigned forwardsClipped = std::distance(adapterRange.adapterRangeEnd_, sequenceEnd);

    clipBackwards = backwardsClipped < forwardsClipped;
    // If we give too much freedom, returned results will provide conflicting alignment positions to
    // the alignment scoring mechanism and it will assign low alignment scores thinking the fragment aligns to a repeat
    const unsigned sequenceLength = std::distance(sequenceBegin, sequenceEnd);
    const std::vector<char> &reference = contig.forward_;
    if (backwardsClipped && forwardsClipped && abs(backwardsClipped - forwardsClipped) < 9)
    {
        // avoid counting matches for alignments that partially fall outside the reference
        if (contigPosition >=0 && reference.size() >= unsigned(contigPosition + sequenceLength))
        {
            std::vector<char>::const_iterator referenceBegin = reference.begin() + contigPosition;
            std::vector<char>::const_iterator referenceEnd = referenceBegin + sequenceLength;
            // count matches if both sides are of close lengths
            const unsigned backwardsMatches = countMatches(
                sequenceBegin, adapterRange.adapterRangeBegin_, referenceBegin, referenceBegin + backwardsClipped);
            const unsigned forwardsMatches = countMatches(
                adapterRange.adapterRangeEnd_, sequenceEnd, referenceEnd - forwardsClipped, referenceEnd);

            ISAAC_THREAD_CERR_DEV_TRACE("backwards matches " << countMismatches(sequenceBegin, adapterRange.adapterRangeBegin_,
                                                                        referenceBegin, referenceBegin + backwardsClipped));
            ISAAC_THREAD_CERR_DEV_TRACE("forwards matches " << countMismatches(adapterRange.adapterRangeEnd_, sequenceEnd,
                                                                        referenceEnd - forwardsClipped, referenceEnd));

            // prefer clipping that leaves more and longer part unclipped
            clipBackwards = (backwardsMatches < forwardsMatches ||
                (backwardsMatches == forwardsMatches && backwardsClipped < forwardsClipped));
        }
    }
    else if (!backwardsClipped || !forwardsClipped)
    {
        // If it is clipping all the way to one of the ends,make sure the side that we clip has decent
        // amount of mismatches. Else, we're clipping something that just happens to have an adapter sequence in it...
        std::vector<char>::const_iterator referenceBegin = reference.begin() + contigPosition;
        std::vector<char>::const_iterator referenceEnd = referenceBegin + sequenceLength;
        if (clipBackwards && !backwardsClipped)
        {
            const unsigned basesClipped = std::distance(sequenceBegin, adapterRange.adapterRangeEnd_);
            ISAAC_THREAD_CERR_DEV_TRACE("backwards percent " << percentMismatches(sequenceBegin, adapterRange.adapterRangeEnd_,
                                                                        referenceBegin, referenceBegin + basesClipped));
            ISAAC_THREAD_CERR_DEV_TRACE("backwards count " << countMismatches(sequenceBegin, adapterRange.adapterRangeEnd_,
                                                                        referenceBegin, referenceBegin + basesClipped));

            return percentMismatches(sequenceBegin, adapterRange.adapterRangeEnd_,
                                     referenceBegin, referenceBegin + basesClipped) > TOO_GOOD_READ_MISMATCH_PERCENT;
        }
        else if (!clipBackwards && !forwardsClipped)
        {
            const unsigned basesClipped = std::distance(adapterRange.adapterRangeBegin_, sequenceEnd);
            ISAAC_THREAD_CERR_DEV_TRACE("forwards percent " << percentMismatches(adapterRange.adapterRangeBegin_, sequenceEnd,
                                                                        referenceEnd - basesClipped, referenceEnd));
            ISAAC_THREAD_CERR_DEV_TRACE("forwards count " << countMismatches(adapterRange.adapterRangeBegin_, sequenceEnd,
                                                                        referenceEnd - basesClipped, referenceEnd));

            return percentMismatches(adapterRange.adapterRangeBegin_, sequenceEnd,
                                     referenceEnd - basesClipped, referenceEnd) > TOO_GOOD_READ_MISMATCH_PERCENT;
        }
    }

    return true;
}
/**
 * \brief Clips off the adapters and part of the sequence that assumedly does not align.
 *        The decision on which part to clip is based on the length.
 *
 * \param contigPosition in/out: on return: new fragmentMetadata.position that points at the
 *                               first not clipped base
 */
void FragmentSequencingAdapterClipper::clip(
    const reference::Contig &contig,
    FragmentMetadata &fragment,
    std::vector<char>::const_iterator &sequenceBegin,
    std::vector<char>::const_iterator &sequenceEnd) const
{
    const SequencingAdapterRange &adapterRange = strandAdapters_.strandRange_[fragment.reverse];
    ISAAC_ASSERT_MSG(adapterRange.initialized_, "checkInitStrand has not been called");
    if (!adapterRange.empty_)
    {
        // non-empty range means adapter sequence identified around the currentBase
        ISAAC_ASSERT_MSG(adapterRange.adapterRangeBegin_ >= sequenceBegin &&
                         adapterRange.adapterRangeBegin_ <= sequenceEnd, "adapter range begin is outside the sequence");
        ISAAC_ASSERT_MSG(adapterRange.adapterRangeEnd_ >= sequenceBegin &&
                         adapterRange.adapterRangeEnd_ <= sequenceEnd, "adapter range end is outside the sequence");

        bool clipBackwards = false;

        if (decideWhichSideToClip(contig, fragment.position, sequenceBegin, sequenceEnd, adapterRange, clipBackwards))
        {
            if (clipBackwards)
            {
                ISAAC_THREAD_CERR_DEV_TRACE("FragmentSequencingAdapterClipper::clip: adapter in " <<
                                            std::string(sequenceBegin, sequenceEnd) <<
                                            ", clipping backwards " <<
                                            std::distance(sequenceBegin, adapterRange.adapterRangeEnd_) << " bases");

                fragment.incrementClipLeft(std::distance(sequenceBegin, adapterRange.adapterRangeEnd_));
                sequenceBegin = adapterRange.adapterRangeEnd_;
            }
            else
            {
                ISAAC_THREAD_CERR_DEV_TRACE("FragmentSequencingAdapterClipper::clip: adapter in " <<
                                            std::string(sequenceBegin, sequenceEnd) <<
                                            ", clipping forwards " <<
                                            std::distance(adapterRange.adapterRangeBegin_, sequenceEnd) << " bases");

                fragment.incrementClipRight(std::distance(adapterRange.adapterRangeBegin_, sequenceEnd));
                sequenceEnd = adapterRange.adapterRangeBegin_;
            }
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE("FragmentSequencingAdapterClipper::clip: adapter in " <<
                                        std::string(sequenceBegin, sequenceEnd) <<
                                        ", not clipping!");
        }
    }
}

} // namespace matchSelector
} // namespace alignment
} // namespace isaac
