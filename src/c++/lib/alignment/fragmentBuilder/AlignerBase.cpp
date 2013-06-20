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
 ** \file AlignerBase.cpp
 **
 ** \brief See AlignerBase.hh
 ** 
 ** \author Roman Petrovski
 **/
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "alignment/fragmentBuilder/AlignerBase.hh"
#include "alignment/Alignment.hh"

namespace isaac
{
namespace alignment
{
namespace fragmentBuilder
{

AlignerBase::AlignerBase(
    const int gapMatchScore,
    const int gapMismatchScore,
    const int gapOpenScore,
    const int gapExtendScore,
    const int minGapExtendScore)
    : normalizedMismatchScore_(gapMatchScore - gapMismatchScore)
    , normalizedGapOpenScore_(gapMatchScore - gapOpenScore)
    , normalizedGapExtendScore_(gapMatchScore - gapExtendScore)
    , normalizedMaxGapExtendScore_(-minGapExtendScore)
{
}

/**
 * \brief Adjusts the sequence iterators to stay within the reference. Adjusts sequenceBeginReferencePosition
 *        to point at the first not clipped base.
 *
 */
void AlignerBase::clipReference(
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
void AlignerBase::clipReadMasking(
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

unsigned AlignerBase::updateFragmentCigar(
    const flowcell::ReadMetadataList &readMetadataList,
    const std::vector<char> &reference,
    FragmentMetadata &fragmentMetadata,
    const long strandPosition,
    const Cigar &cigarBuffer,
    const unsigned cigarOffset) const
{
    const Read &read = fragmentMetadata.getRead();
    const bool reverse = fragmentMetadata.reverse;
    const std::vector<char> &sequence = read.getStrandSequence(reverse);
    const std::vector<char> &quality = read.getStrandQuality(reverse);

    ISAAC_ASSERT_MSG(!reference.empty(), "Reference contig was not loaded for " << fragmentMetadata);

    ISAAC_ASSERT_MSG(0 <= strandPosition, "position must be positive for CIGAR update " << fragmentMetadata << " strandPosition:" << strandPosition <<
                     " CIGAR: " << Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarBuffer.end()));
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
            unsigned matchesInARow = 0;
            for (unsigned j = 0; length > j; ++j)
            {
                if (isMatch(sequence[currentBase], *currentReference))
                {
                    ++matchCount;
                    ++matchesInARow;
                    fragmentMetadata.logProbability += Quality::getLogMatch(quality[currentBase]);
                }
                else
                {
                    fragmentMetadata.matchesInARow = std::max(fragmentMetadata.matchesInARow, matchesInARow);
                    matchesInARow = 0;
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
            fragmentMetadata.matchesInARow = std::max(fragmentMetadata.matchesInARow, matchesInARow);
        }
        else if (opCode == Cigar::INSERT)
        {
            currentBase += length;
            fragmentMetadata.editDistance += length;
            ++fragmentMetadata.gapCount;
            fragmentMetadata.smithWatermanScore += normalizedGapOpenScore_ + std::min(normalizedMaxGapExtendScore_, (length - 1) * normalizedGapExtendScore_);
        }
        else if (opCode == Cigar::DELETE)
        {
            currentReference += length;
            fragmentMetadata.editDistance += length;
            ++fragmentMetadata.gapCount;
            fragmentMetadata.smithWatermanScore += normalizedGapOpenScore_ + std::min(normalizedMaxGapExtendScore_, (length - 1) * normalizedGapExtendScore_);
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
    fragmentMetadata.observedLength = currentReference - reference.begin() - strandPosition;
    fragmentMetadata.position = strandPosition;
    ISAAC_ASSERT_MSG(currentBase == sequence.size(),
                     "Unexpected discrepancy between cigar and sequence" << fragmentMetadata);

    return matchCount;
}

} // namespace fragmentBuilder
} // namespace alignment
} // namespace isaac
