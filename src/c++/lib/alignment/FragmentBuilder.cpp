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
 ** \file FragmentBuilder.cpp
 **
 ** \brief See FragmentBuilder.hh
 ** 
 ** \author Come Raczy
 **/

#include <numeric>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "alignment/Alignment.hh"
#include "alignment/FragmentBuilder.hh"
#include "alignment/Quality.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

FragmentBuilder::FragmentBuilder(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const unsigned repeatThreshold,
    const unsigned maxSeedsPerRead,
    const unsigned gappedMismatchesMax,
    const int gapMatchScore,
    const int gapMismatchScore,
    const int gapOpenScore,
    const int gapExtendScore)
    : normalizedMismatchScore_(gapMatchScore - gapMismatchScore)
    , normalizedGapOpenScore_(gapMatchScore - gapOpenScore)
    , normalizedGapExtendScore_(gapMatchScore - gapExtendScore)
    , repeatThreshold_(repeatThreshold)
    , gappedMismatchesMax_(gappedMismatchesMax)
    , bandedSmithWaterman_(gapMatchScore, gapMismatchScore, -gapOpenScore, -gapExtendScore,
                           flowcell::getMaxTotalReadLength(flowcellLayoutList))
    , seedMatchCounts_(maxSeedsPerRead * readsMax_)
    , repeatSeedsCount_(0)
    , fragments_(readsMax_) // max number of read ends we ever have to deal with
    , cigarBuffer_(Cigar::getMaxOperationsForReads(flowcellLayoutList) *
                   // one seed generates up to repeat threshold matches for each strand
                   repeatThreshold_ * maxSeedsPerRead * 2)
{
    std::for_each(fragments_.begin(), fragments_.end(),
                  boost::bind(&std::vector<FragmentMetadata>::reserve, _1,
                              // one seed generates up to repeat threshold matches for each strand
                              repeatThreshold_ * maxSeedsPerRead * 2));
}

void FragmentBuilder::clear()
{
    BOOST_FOREACH(std::vector<FragmentMetadata> &fragmentList, fragments_)
    {
        fragmentList.clear();
    }
    cigarBuffer_.clear();
    std::fill(seedMatchCounts_.begin(), seedMatchCounts_.end(), 0);
    repeatSeedsCount_ = 0;
}

/**
 * \return true if at least one fragment was built.
 */
bool FragmentBuilder::build(
    const std::vector<reference::Contig> &contigList,
    const flowcell::ReadMetadataList &readMetadataList,
    const SeedMetadataList &seedMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    std::vector<Match>::const_iterator matchBegin,
    const std::vector<Match>::const_iterator matchEnd,
    const Cluster &cluster,
    const bool withGaps)
{
    clear();
    if (matchBegin < matchEnd)
    {
        ISAAC_ASSERT_MSG(!matchBegin->isNoMatch(), "Fake match lists must be dealt with outside");
        ISAAC_ASSERT_MSG(cluster.getNonEmptyReadsCount() == readMetadataList.size(), "cluster geometry must match");
        ISAAC_THREAD_CERR_DEV_TRACE("FragmentBuilder::build: cluster " << matchBegin->seedId.getCluster() << " (" << cluster.getId() << ")");

        for(;matchEnd != matchBegin && !matchBegin->isNoMatch(); ++matchBegin)
        {
            if (repeatThreshold_ > seedMatchCounts_[matchBegin->getSeedId().getSeed()])
            {
                if (matchBegin->isTooManyMatch())
                {
                    seedMatchCounts_[matchBegin->getSeedId().getSeed()] = repeatThreshold_;
                    ++repeatSeedsCount_;
                    const unsigned matchReadIndex = seedMetadataList[matchBegin->getSeedId().getSeed()].getReadIndex();
                    ISAAC_ASSERT_MSG(fragments_[matchReadIndex].empty(), "Too-many matches are expected to sort to the top");
                    ISAAC_THREAD_CERR_DEV_TRACE("FragmentBuilder::build: blocking read " << matchReadIndex  << ", seed " << matchBegin->getSeedId().getSeed());

                }
                else
                {
                    if (repeatThreshold_ == ++seedMatchCounts_[matchBegin->getSeedId().getSeed()])
                    {
                        ++repeatSeedsCount_;
                        ISAAC_THREAD_CERR_DEV_TRACE("FragmentBuilder::build: read " << seedMetadataList[matchBegin->getSeedId().getSeed()].getReadIndex()  <<
                                                    ", seed " << matchBegin->getSeedId().getSeed() << " exceeded repeat threshold and will be ignored");
                    }
                    else
                    {
                        addMatch(readMetadataList, seedMetadataList, *matchBegin, cluster);
                    }
                }
            }
        }

        if (repeatSeedsCount_)
        {
            // if some of the seeds produced more than repeatThreshold_ matches,
            // we need to remove those matches to avoid bias
            std::for_each(fragments_.begin(), fragments_.end(),
                          boost::bind(&FragmentBuilder::removeRepeatSeedAlignments, this, _1));
        }

        if (fragments_.end() != std::find_if(fragments_.begin(), fragments_.end(),
                                             !boost::bind(&std::vector<FragmentMetadata>::empty, _1)))
        {
            std::for_each(fragments_.begin(), fragments_.end(), consolidateDuplicateFragments);
            alignFragments(contigList, readMetadataList, sequencingAdapters, withGaps);

            // gapped alignment and adapter trimming may adjust the alignment position
            std::for_each(fragments_.begin(), fragments_.end(), consolidateDuplicateFragments);
            return true;
        }
    }
    return false;
}

void FragmentBuilder::alignFragments(
    const std::vector<reference::Contig> &contigList,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const bool withGaps)
{
    unsigned fragmentIndex = 0;
    BOOST_FOREACH(std::vector<FragmentMetadata> &fragmentList, fragments_)
    {
        if (!fragmentList.empty())
        {
            ISAAC_DEV_TRACE_BLOCK(const Cluster & cluster = fragmentList.at(0).getCluster();)
            ISAAC_THREAD_CERR_DEV_TRACE("Listing fragments for cluster " << cluster.getTile() <<
                ":" << cluster.getId() << " [" << fragmentIndex << "] " << sequencingAdapters.size() << " adapters");
                
            matchSelector::FragmentSequencingAdapterClipper adapterClipper(sequencingAdapters);
            BOOST_FOREACH(FragmentMetadata &fragmentMetadata, fragmentList)
            {
                fragmentMetadata.repeatSeedsCount = repeatSeedsCount_;
                ISAAC_THREAD_CERR_DEV_TRACE("    Original   : " << fragmentMetadata);
                const unsigned contigId = fragmentMetadata.contigId;
                ISAAC_ASSERT_MSG(contigList.size() > contigId, "Unexpected contig id");
                const reference::Contig &contig = contigList[contigId];

                adapterClipper.checkInitStrand(fragmentMetadata, contig);
                alignUngapped(fragmentMetadata, cigarBuffer_, readMetadataList, adapterClipper, contig);
                ISAAC_THREAD_CERR_DEV_TRACE("    Aligned    : " << fragmentMetadata);
                if (withGaps && BandedSmithWaterman::mismatchesCutoff < fragmentMetadata.mismatchCount)
                {
                    FragmentMetadata tmp = fragmentMetadata;
                    const unsigned matchCount = alignGapped(tmp, cigarBuffer_, readMetadataList, adapterClipper, contig);
                    ISAAC_THREAD_CERR_DEV_TRACE("    Gap-aligned: " << tmp);
                    if (matchCount && matchCount + bandedSmithWaterman_.widestGapSize > fragmentMetadata.getObservedLength() &&
                        (tmp.mismatchCount <= gappedMismatchesMax_) &&
                        (fragmentMetadata.mismatchCount > tmp.mismatchCount) &&
                        ISAAC_LP_LESS(fragmentMetadata.logProbability, tmp.logProbability))
                    {
                        ISAAC_THREAD_CERR_DEV_TRACE("    Using gap-aligned: " << tmp);
                        fragmentMetadata = tmp;
                    }
                }
            }
        }
        ++fragmentIndex;
    }
}

void FragmentBuilder::addMatch(
    const flowcell::ReadMetadataList &readMetadataList,
    const SeedMetadataList &seedMetadataList, const Match &match,
    const Cluster &cluster)
{
    const unsigned seedIndex = match.seedId.getSeed();
    const SeedMetadata &seedMetadata = seedMetadataList.at(seedIndex);
    ISAAC_ASSERT_MSG(seedMetadata.getLength() == CURRENTLY_SUPPORTED_SEED_LENGTH,
                     "Currently only 32-base seeds are supported.");
    const unsigned readIndex = seedMetadata.getReadIndex();
    const reference::ReferencePosition &seedLocation = match.location;
    const unsigned contigId = seedLocation.getContigId();
    const bool reverse = match.seedId.isReverse();
    const long readPosition = getReadPosition(readMetadataList, seedMetadata, seedLocation.getPosition(), reverse);

    ISAAC_THREAD_CERR_DEV_TRACE("    adding match at position " << contigId << ":" << readPosition);
    // Cigar buffer and offset are set separately at the point when we align
    fragments_[readIndex].push_back(FragmentMetadata(&cluster, NULL, readIndex));
    FragmentMetadata &fragment = fragments_[readIndex].back();
    fragment.firstSeedIndex = seedIndex;
    fragment.contigId = contigId;
    fragment.position = readPosition;
    fragment.reverse = reverse;
    if (seedLocation.hasNeighbors())
    {
        fragment.nonUniqueSeedOffsets.push_back(seedMetadata.getOffset());
    }
    else
    {
        fragment.uniqueSeedCount = 1;
    }
}


inline bool exceeedsThreshold(const std::vector<unsigned> &vect, const unsigned index, const unsigned threshold)
{
    return vect[index] >= threshold;
}
/**
 * \brief Removes entries introduced by the seeds which later on have been found to have over repeatThreshold_
 *        matches
 *
 */
void FragmentBuilder::removeRepeatSeedAlignments(std::vector<FragmentMetadata> &fragmentList)
{
    fragmentList.erase(std::remove_if(fragmentList.begin(), fragmentList.end(),
                                      boost::bind(&exceeedsThreshold, boost::ref(seedMatchCounts_), boost::bind(&FragmentMetadata::getFirstSeedIndex, _1), repeatThreshold_)),
                                      fragmentList.end());
}

/**
 * \brief Removes entries designating the same alignment location and strand
 *
 * This is important not only to avoid repetitive processing but also to avoid good alignments
 * scored poorly because they are present multiple times in the list.
 */
void FragmentBuilder::consolidateDuplicateFragments(std::vector<FragmentMetadata> &fragmentList)
{
    if (fragmentList.size() < 2)
    {
        return;
    }
    // although initially matches arrive ordered by location, gapped alignment might have moved the
    // start position of some.
    std::sort(fragmentList.begin(), fragmentList.end());
    std::vector<FragmentMetadata>::iterator lastFragment = fragmentList.begin();
    using boost::format;
    ISAAC_THREAD_CERR_DEV_TRACE((format("    adding fragment at %d:%d") %  lastFragment->contigId % lastFragment->position).str());
    for (std::vector<FragmentMetadata>::iterator currentFragment = lastFragment + 1; fragmentList.end() != currentFragment; ++currentFragment)
    {
        if (lastFragment->contigId == currentFragment->contigId &&
            lastFragment->position == currentFragment->position &&
            lastFragment->reverse == currentFragment->reverse)
        {
            ISAAC_THREAD_CERR_DEV_TRACE((format("    adding %2d unique seeds to %d:%d") %  currentFragment->uniqueSeedCount %
                                         currentFragment->contigId % currentFragment->position).str());
            lastFragment->uniqueSeedCount += currentFragment->uniqueSeedCount;
            std::copy(currentFragment->nonUniqueSeedOffsets.begin(), currentFragment->nonUniqueSeedOffsets.end(),
                      std::back_inserter(lastFragment->nonUniqueSeedOffsets));
        }
        else
        {
            ++lastFragment;
            ISAAC_THREAD_CERR_DEV_TRACE((format("    adding fragment at %d:%d") %  lastFragment->contigId % lastFragment->position).str());
            if (lastFragment != currentFragment)
            {
                *lastFragment = *currentFragment;
            }
        }
    }
    fragmentList.resize(1 + lastFragment - fragmentList.begin());
}

inline long FragmentBuilder::getReadPosition(
    const flowcell::ReadMetadataList &readMetadataList,
    const SeedMetadata &seedMetadata, const long seedPosition, const bool reverse) const
{
    const int seedOffset = seedMetadata.getOffset();
    if (reverse)
    {
        const unsigned readIndex = seedMetadata.getReadIndex();
        const unsigned readLength = readMetadataList.at(readIndex).getLength();
        const unsigned seedLength = seedMetadata.getLength();
        // 'seedPosition + seedLength + seedOffset' is the first position past the end of the read
        return seedPosition + seedLength + seedOffset - readLength;
    }
    else
    {
        return seedPosition - seedOffset;
    }
}

/**
 * \brief Adjusts the sequence iterators to stay within the reference. Adjusts sequenceBeginReferencePosition
 *        to point at the first not clipped base.
 *
 */
static void clipReference(
    const long referenceSize,
    FragmentMetadata &fragment,
    std::vector<char>::const_iterator &sequenceBegin,
    std::vector<char>::const_iterator &sequenceEnd)
{
    const long referenceLeft = referenceSize - fragment.position;
    if (referenceLeft >= 0)
    {
        if (referenceLeft < std::distance(sequenceBegin, sequenceEnd))
        {
            sequenceEnd = sequenceBegin + referenceLeft;
        }

        if (0 > fragment.position)
        {
            sequenceBegin -= fragment.position;
            fragment.position = 0L;
        }

        // in some cases other clipping can end the sequence before the reference even begins
        // or begin after it ends...
        sequenceEnd = std::max(sequenceEnd, sequenceBegin);
    }
    else
    {
        // the picard sam ValidateSamFile does not like it when alignment position points to the next base after the end of the contig.
        fragment.position += referenceLeft - 1;
        sequenceBegin += referenceLeft - 1;
        --sequenceBegin;
        sequenceEnd = sequenceBegin;
    }
}

/**
 * \brief Sets the sequence iterators according to the masking information stored in the read.
 *        Adjusts fragment.position to point at the first non-clipped base.
 *
 */
static void clipReadMasking(
    const alignment::Read &read,
    FragmentMetadata &fragment,
    std::vector<char>::const_iterator &sequenceBegin,
    std::vector<char>::const_iterator &sequenceEnd)
{
    std::vector<char>::const_iterator maskedBegin;
    std::vector<char>::const_iterator maskedEnd;
    if (fragment.reverse)
    {
        maskedBegin = read.getReverseSequence().begin() + read.getEndCyclesMasked();
        maskedEnd = read.getReverseSequence().end() - read.getBeginCyclesMasked();
    }
    else
    {
        maskedBegin = read.getForwardSequence().begin() + read.getBeginCyclesMasked();
        maskedEnd = read.getForwardSequence().end() - read.getEndCyclesMasked();
    }

    if (maskedBegin > sequenceBegin)
    {
        fragment.incrementClipLeft(std::distance(sequenceBegin, maskedBegin));
        sequenceBegin = maskedBegin;
    }

    if (maskedEnd < sequenceEnd)
    {
        fragment.incrementClipRight(std::distance(maskedEnd, sequenceEnd));
        sequenceEnd = maskedEnd;
    }
}

void FragmentBuilder::alignUngapped(
    FragmentMetadata &fragmentMetadata,
    Cigar &cigarBuffer,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::FragmentSequencingAdapterClipper &adapterClipper,
    const reference::Contig &contig) const
{
    const unsigned cigarOffset = cigarBuffer.size();

    fragmentMetadata.resetAlignment(cigarBuffer);

    const std::vector<char> &reference = contig.forward_;

    const Read &read = fragmentMetadata.getRead();

    const bool reverse = fragmentMetadata.reverse;
    std::vector<char>::const_iterator sequenceBegin = read.getStrandSequence(reverse).begin();
    std::vector<char>::const_iterator sequenceEnd = read.getStrandSequence(reverse).end();

    adapterClipper.clip(contig, fragmentMetadata, sequenceBegin, sequenceEnd);
    clipReadMasking(read, fragmentMetadata, sequenceBegin, sequenceEnd);

    clipReference(reference.size(), fragmentMetadata, sequenceBegin, sequenceEnd);

    const std::vector<char> &sequence = read.getStrandSequence(reverse);

    const unsigned firstMappedBaseOffset = std::distance(sequence.begin(), sequenceBegin);
    if (firstMappedBaseOffset)
    {
        cigarBuffer.addOperation(firstMappedBaseOffset, Cigar::SOFT_CLIP);
    }

    const unsigned mappedBases = std::distance(sequenceBegin, sequenceEnd);
    if (mappedBases)
    {
        const Cigar::OpCode opCode = Cigar::ALIGN;
        cigarBuffer.addOperation(mappedBases, opCode);
    }

    const unsigned clipEndBases = std::distance(sequenceEnd, sequence.end());
    if (clipEndBases)
    {
        cigarBuffer.addOperation(clipEndBases, Cigar::SOFT_CLIP);
    }

    updateFragmentCigar(readMetadataList, reference, fragmentMetadata,
                        fragmentMetadata.position, cigarBuffer, cigarOffset);
}

/// calculate the left and right flanks of the database WRT the query
static std::pair<unsigned, unsigned> getFlanks(
    const long strandPosition,
    const unsigned readLength,
    const unsigned long referenceSize,
    const unsigned widestGapSize)
{
    assert(widestGapSize >= 2);
    if (strandPosition >= widestGapSize/2)
    {
        // take into account the rounding
        if (strandPosition + readLength + (widestGapSize - widestGapSize/2) < static_cast<long>(referenceSize))
        {
            const unsigned left = widestGapSize/2;
            const unsigned right = widestGapSize - left - 1;
            return std::pair<unsigned, unsigned>(left, right);
        }
        else
        {
            assert((long)referenceSize >= readLength + strandPosition);
            const unsigned right = referenceSize - readLength - strandPosition;
            const unsigned left = widestGapSize - right - 1;
            return std::pair<unsigned, unsigned>(left, right);
        }
    }
    else
    {
        assert(strandPosition >= 0);
        const unsigned left = strandPosition;
        const unsigned right = widestGapSize - left - 1;
        return std::pair<unsigned, unsigned>(left, right);
    }
}

unsigned FragmentBuilder::updateFragmentCigar(
    const flowcell::ReadMetadataList &readMetadataList,
    const std::vector<char> &reference,
    FragmentMetadata &fragmentMetadata,
    long strandPosition,
    Cigar &cigarBuffer,
    const unsigned cigarOffset) const
{
    const Read &read = fragmentMetadata.getRead();
    const bool reverse = fragmentMetadata.reverse;
    const std::vector<char> &sequence = read.getStrandSequence(reverse);
    const std::vector<char> &quality = read.getStrandQuality(reverse);

    assert(0 <= strandPosition);
    std::vector<char>::const_iterator currentReference = reference.begin() + strandPosition;

    const unsigned firstCycle = readMetadataList[fragmentMetadata.readIndex].getFirstCycle();
    const unsigned lastCycle = readMetadataList[fragmentMetadata.readIndex].getLastCycle();

    fragmentMetadata.cigarBuffer = &cigarBuffer;
    fragmentMetadata.cigarOffset = cigarOffset;
    fragmentMetadata.cigarLength = cigarBuffer.size() - fragmentMetadata.cigarOffset;
    // adjust cigarOffset and cigarLength
    ISAAC_ASSERT_MSG(cigarBuffer.size() > fragmentMetadata.cigarOffset, "Expecting the new cigar is not empty");

    unsigned currentBase = 0;
    unsigned matchCount = 0;
    for (unsigned i = 0; fragmentMetadata.cigarLength > i; ++i)
    {
        const std::pair<unsigned, Cigar::OpCode> cigar = Cigar::decode(cigarBuffer[fragmentMetadata.cigarOffset + i]);
        const unsigned length = cigar.first;
        const Cigar::OpCode opCode = cigar.second;
        if (opCode == Cigar::ALIGN)
        {
            for (unsigned j = 0; length > j; ++j)
            {
                if (isMatch(sequence[currentBase], *currentReference))
                {
                    ++matchCount;
                    fragmentMetadata.logProbability += Quality::getLogMatch(quality[currentBase]);
                }
                else
                {
                    fragmentMetadata.addMismatchCycle(reverse ? lastCycle - currentBase : firstCycle + currentBase);
                    fragmentMetadata.logProbability += Quality::getLogMismatchFast(quality[currentBase]);
                    fragmentMetadata.smithWatermanScore += normalizedMismatchScore_;
                }
                // the edit distance includes all mismatches and ambiguous bases (Ns)
                if (sequence[currentBase] != *currentReference)
                {
                    ++fragmentMetadata.editDistance;
                }
                ++currentReference;
                ++currentBase;
            }
        }
        else if (opCode == Cigar::INSERT)
        {
            currentBase += length;
            fragmentMetadata.editDistance += length;
            ++fragmentMetadata.gapCount;
            fragmentMetadata.smithWatermanScore += normalizedGapOpenScore_ + (length - 1) * normalizedGapExtendScore_;
        }
        else if (opCode == Cigar::DELETE)
        {
            currentReference += length;
            fragmentMetadata.editDistance += length;
            ++fragmentMetadata.gapCount;
            fragmentMetadata.smithWatermanScore += normalizedGapOpenScore_ + (length - 1) * normalizedGapExtendScore_;
        }
        else if (opCode == Cigar::SOFT_CLIP)
        {
            ISAAC_ASSERT_MSG(0 == i || i + 1 == fragmentMetadata.cigarLength, "Soft clippings are expected to be "
                "found only at the ends of cigar string");
            using boost::lambda::bind;
            using boost::lambda::_1;
            using boost::lambda::_2;
            fragmentMetadata.logProbability =
                std::accumulate(quality.begin() + currentBase, quality.begin() + currentBase + length,
                                fragmentMetadata.logProbability,
                                bind(std::plus<double>(), _1, bind(Quality::getLogMatch, _2)));

            // NOTE! Not advancing the reference for soft clips
            currentBase += length;
        }
        else
        {
            using boost::format;
            const format message = format("Unexpected Cigar OpCode: %d") % opCode;
            BOOST_THROW_EXCEPTION(common::PostConditionException(message.str()));
        }
    }
    ISAAC_ASSERT_MSG(currentBase == sequence.size(), "Unexpected discrepancy between cigar and sequence");
    fragmentMetadata.observedLength = currentReference - reference.begin() - strandPosition;
    fragmentMetadata.position = strandPosition;

    return matchCount;
}

unsigned FragmentBuilder::alignGapped(
    FragmentMetadata &fragmentMetadata,
    Cigar &cigarBuffer,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::FragmentSequencingAdapterClipper &adapterClipper,
    const reference::Contig &contig) const
{
    const unsigned cigarOffset = cigarBuffer.size();
    fragmentMetadata.resetAlignment(cigarBuffer);

    const Read &read = fragmentMetadata.getRead();
    const bool reverse = fragmentMetadata.reverse;
    const std::vector<char> &sequence = read.getStrandSequence(reverse);
    const std::vector<char> &reference = contig.forward_;

    std::vector<char>::const_iterator sequenceBegin = sequence.begin();
    std::vector<char>::const_iterator sequenceEnd = sequence.end();

    adapterClipper.clip(contig, fragmentMetadata, sequenceBegin, sequenceEnd);
    clipReadMasking(read, fragmentMetadata, sequenceBegin, sequenceEnd);

    clipReference(reference.size(), fragmentMetadata, sequenceBegin, sequenceEnd);

    const unsigned firstMappedBaseOffset = std::distance(sequence.begin(), sequenceBegin);
    if (firstMappedBaseOffset)
    {
        cigarBuffer.addOperation(firstMappedBaseOffset, Cigar::SOFT_CLIP);
    }

    const unsigned sequenceLength = std::distance(sequenceBegin, sequenceEnd);

    // position of the fragment on the strand
    long strandPosition = fragmentMetadata.position;
    ISAAC_ASSERT_MSG(0 <= strandPosition, "alignUngapped should have clipped reads beginning before the reference");
    // no gapped alignment if the reference is too short
    const unsigned widestGapSize = bandedSmithWaterman_.widestGapSize;

    if (static_cast<long>(reference.size()) < sequenceLength + strandPosition + widestGapSize)
    {
        ISAAC_THREAD_CERR_DEV_TRACE("alignGapped: reference too short!");
        return 0;
    }
    // find appropriate beginning and end for the database
    const std::pair<unsigned, unsigned> flanks = getFlanks(strandPosition, sequenceLength, reference.size(), widestGapSize);
    assert(flanks.first + flanks.second == widestGapSize - 1);
    assert(flanks.first <= strandPosition);
    assert(strandPosition + sequenceLength + flanks.second <= (long)reference.size());
    const std::vector<char>::const_iterator databaseBegin = reference.begin() + strandPosition - flanks.first;
    const std::vector<char>::const_iterator databaseEnd = databaseBegin + flanks.first + sequenceLength + flanks.second;

    ISAAC_THREAD_CERR_DEV_TRACE("Gap-aligning " << std::string(sequenceBegin, sequenceEnd) <<
        " against " << std::string(databaseBegin, databaseEnd));
    strandPosition += bandedSmithWaterman_.align(sequenceBegin, sequenceEnd, databaseBegin, databaseEnd, cigarBuffer);

    const unsigned clipEndBases = std::distance(sequenceEnd, sequence.end());
    if (clipEndBases)
    {
        cigarBuffer.addOperation(clipEndBases, Cigar::SOFT_CLIP);
    }

    // adjust the start position of the fragment
    strandPosition -= flanks.first;

    const unsigned matchCount = updateFragmentCigar(readMetadataList, reference, fragmentMetadata,
                                                    strandPosition, cigarBuffer, cigarOffset);

    return matchCount;
}

} // namespace alignment
} // namespace isaac
