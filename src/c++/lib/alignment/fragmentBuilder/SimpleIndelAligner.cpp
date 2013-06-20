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
 ** \file SimpleIndelAligner.cpp
 **
 ** \brief See SimpleIndelAligner.hh
 ** 
 ** \author Roman Petrovski
 **/

#include "alignment/fragmentBuilder/SimpleIndelAligner.hh"
#include "alignment/Alignment.hh"

namespace isaac
{
namespace alignment
{
namespace fragmentBuilder
{
const unsigned SimpleIndelAligner::GAP_FLANK_BASES;

SimpleIndelAligner::SimpleIndelAligner(
    const int gapMatchScore,
    const int gapMismatchScore,
    const int gapOpenScore,
    const int gapExtendScore,
    const int minGapExtendScore,
    const unsigned semialignedGapLimit)
    : AlignerBase(gapMatchScore, gapMismatchScore, gapOpenScore, gapExtendScore, minGapExtendScore)
    , semialignedGapLimit_(semialignedGapLimit)
{
}

/**
 * \brief Patches the front fragment with cigar that produces the lowest number of mismatches assuming there
 *        is a deletion in the read somewhere between the frontFragment first seed and back fragment first seed
 *
 * \param headSeedOffset offset of the seed for headAlignment from the left-most base in respect to the reference
 * \param tailSeedOffset offset of the seed for tailAlignment from the left-most base in respect to the reference
 */
void SimpleIndelAligner::alignSimpleDeletion(
    Cigar &cigarBuffer,
    FragmentMetadata &headAlignment,
    const unsigned headSeedOffset,
    FragmentMetadata &tailAlignment,
    const unsigned tailSeedOffset,
    const unsigned tailSeedLength,
    const std::vector<reference::Contig> &contigList,
    const flowcell::ReadMetadataList &readMetadataList) const
{
    if (headSeedOffset < headAlignment.getBeginClippedLength())
    {
        ISAAC_THREAD_CERR_DEV_TRACE("alignSimpleDeletion: clipping on the left overlaps the head seed.");
        return;
    }

    if (tailAlignment.getBeginClippedLength() + tailAlignment.getObservedLength() < tailSeedOffset + tailSeedLength)
    {
        ISAAC_THREAD_CERR_DEV_TRACE("alignSimpleDeletion: clipping on the right overlaps the tail seed. "
            "tailSeedOffset + tailSeedLength:" << tailSeedOffset  + tailSeedLength <<
            " tailAlignment.getBeginClippedLength() + tailAlignment.getObservedLength():" << tailAlignment.getBeginClippedLength() + tailAlignment.getObservedLength());
        return;
    }

    ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleDeletion: " << headAlignment << "-" << tailAlignment);

    // allow deletion to be placed within the first seed to capture the earliest possible location.
    const unsigned tailOffset = headSeedOffset;

    const Read &read = headAlignment.getRead();
    const bool reverse = headAlignment.reverse;
    const std::vector<char>::const_iterator sequenceBegin = read.getStrandSequence(reverse).begin();
    const std::vector<char> &reference = contigList[headAlignment.contigId].forward_;

    std::vector<char>::const_iterator tailIterator = sequenceBegin + tailOffset;
    unsigned tailLength = tailAlignment.getBeginClippedLength() + tailAlignment.getObservedLength() - tailOffset;
    // number of tail mismatches when deletion is not present
    const unsigned tailMismatches = countMismatches(tailIterator,
                                                    reference.begin() + headAlignment.getUnclippedPosition() + tailOffset, reference.end(),
                                                    tailLength, &boost::cref<char>);
    if (!tailMismatches)
    {
        ISAAC_THREAD_CERR_DEV_TRACE("alignSimpleDeletion: no point to try, the head alignment is already good enough");
        return;
    }

    const unsigned deletionLength = boost::numeric_cast<unsigned>(tailAlignment.getUnclippedPosition() - headAlignment.getUnclippedPosition());

    // try to introduce the deletion on the headFragment see if the number of mismatches reduces below the original

    // number of mismatches when deletion is at the leftmost possible position
    unsigned rightRealignedMismatches = countMismatches(tailIterator,
                                                        reference.begin() + tailAlignment.getUnclippedPosition() + tailOffset, reference.end(),
                                                        tailLength, &boost::cref<char>);
    // we're starting at the situation where the whole tail of the head alignment is moved by deletionLength
    unsigned leftRealignedMismatches = 0;
    unsigned leftFlankMismatches = countMismatches(tailIterator - std::min(GAP_FLANK_BASES, tailOffset),
                                                   reference.begin() + headAlignment.getUnclippedPosition() + tailOffset - std::min(32U, tailOffset), reference.end(),
                                                   std::min(GAP_FLANK_BASES, tailOffset), &boost::cref<char>);;

    unsigned rightFlankMismatches = countMismatches(tailIterator,
                                                    reference.begin() + tailAlignment.getUnclippedPosition() + tailOffset, reference.end(),
                                                    std::min(GAP_FLANK_BASES, tailLength), &boost::cref<char>);;

    ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleDeletion " <<
                                tailMismatches << "htmm " << rightRealignedMismatches << ":" << leftRealignedMismatches << "rhtrmm:lhtrmm ");

    std::vector<char>::const_iterator referenceIterator = reference.begin() + headAlignment.getUnclippedPosition() + tailOffset;

    unsigned bestMismatches = tailMismatches;
    unsigned bestLeftFlankMismatches = leftFlankMismatches;
    unsigned bestRightFlankMismatches = rightFlankMismatches;
    unsigned bestOffset = -1U;

    ISAAC_THREAD_CERR_DEV_TRACE("headTailOffset=" << tailOffset << " tailSeedOffset=" << tailSeedOffset << " deletionLength=" << deletionLength);

    for (unsigned deletionOffset = tailOffset; bestMismatches && deletionOffset <= tailSeedOffset;
        ++deletionOffset, ++tailIterator, ++referenceIterator, --tailLength)
    {
        const unsigned thisOffsetMismatches = leftRealignedMismatches + rightRealignedMismatches;

        if (bestMismatches > thisOffsetMismatches)
        {
            bestOffset = deletionOffset;
            bestMismatches = thisOffsetMismatches;
            bestLeftFlankMismatches = leftFlankMismatches;
            bestRightFlankMismatches = rightFlankMismatches;
        }

        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(headAlignment.getCluster().getId(), "bestMismatches=" << bestMismatches <<
                                    " leftRealignedMismatches=" << leftRealignedMismatches <<
                                    " rightRealignedMismatches=" << rightRealignedMismatches << " deletionOffset=" << deletionOffset <<
                                    " leftFlankMismatches=" << leftFlankMismatches <<
                                    " rightFlankMismatches=" << rightFlankMismatches);
//        ISAAC_THREAD_CERR_DEV_TRACE(" headTail=" << std::string(tailIterator, read.getStrandSequence(reverse).end()) <<
//                                    " headTailRef=" << std::string(referenceIterator + deletionLength, referenceIterator + deletionLength + std::distance(tailIterator, read.getStrandSequence(reverse).end())));
        ISAAC_ASSERT_MSG(rightRealignedMismatches != -1U, "tada");
        const bool newLeftMismatch = !isMatch(*tailIterator, *referenceIterator);
        leftRealignedMismatches += newLeftMismatch;
        leftFlankMismatches += newLeftMismatch;
        if (deletionOffset >= GAP_FLANK_BASES)
        {
            leftFlankMismatches -= !isMatch(*(tailIterator - GAP_FLANK_BASES), *(referenceIterator - GAP_FLANK_BASES));
        }
        const bool disappearingRightMismatch = !isMatch(*tailIterator, *(referenceIterator + deletionLength));
        rightRealignedMismatches -= disappearingRightMismatch;
        rightFlankMismatches -= disappearingRightMismatch;
        if (tailLength > GAP_FLANK_BASES)
        {
            rightFlankMismatches += !isMatch(*(tailIterator + GAP_FLANK_BASES), *(referenceIterator + deletionLength + GAP_FLANK_BASES));
        }
    }


    if (bestLeftFlankMismatches <= GAP_FLANK_MISMATCHES_MAX && bestRightFlankMismatches <= GAP_FLANK_MISMATCHES_MAX && -1U != bestOffset)
    {
        const long clippingPositionOffset = headAlignment.getBeginClippedLength();
        const unsigned leftMapped = bestOffset - clippingPositionOffset;

        const unsigned headMismatches = countMismatches(sequenceBegin + clippingPositionOffset,
                                                        reference.begin() + headAlignment.position, reference.end(),
                                                        leftMapped, &boost::cref<char>);

        const unsigned newMismatches = headMismatches + bestMismatches;

        const unsigned sws = normalizedMismatchScore_ * newMismatches + normalizedGapOpenScore_ +
            std::min(normalizedMaxGapExtendScore_, (deletionLength - 1) * normalizedGapExtendScore_);
        if (headAlignment.smithWatermanScore > sws || (headAlignment.smithWatermanScore == sws && headAlignment.getMismatchCount() > newMismatches))
        {
            const unsigned cigarOffset = cigarBuffer.size();

            if (clippingPositionOffset)
            {
                cigarBuffer.addOperation(clippingPositionOffset, Cigar::SOFT_CLIP);
            }

            if (leftMapped)
            {
                cigarBuffer.addOperation(leftMapped, Cigar::ALIGN);
                cigarBuffer.addOperation(deletionLength, Cigar::DELETE);
            }
            else
            {
                // prevent cigars starting from deletion.
                headAlignment.position += deletionLength;
            }

            const unsigned rightMapped = headAlignment.getObservedLength()  + headAlignment.getEndClippedLength() - leftMapped - tailAlignment.getEndClippedLength();
            if (rightMapped)
            {
                cigarBuffer.addOperation(rightMapped, Cigar::ALIGN);
            }

            const unsigned clipEndBases = tailAlignment.getEndClippedLength();
            if (clipEndBases)
            {
                cigarBuffer.addOperation(clipEndBases, Cigar::SOFT_CLIP);
            }

            headAlignment.resetAlignment(cigarBuffer);
            // carry over alignment-independent clipping information (quality trimming, adapter masking)
            headAlignment.rightClipped() = tailAlignment.rightClipped();
            ISAAC_ASSERT_MSG(updateFragmentCigar(readMetadataList, reference, headAlignment,
                                                 headAlignment.position + clippingPositionOffset, cigarBuffer, cigarOffset),
                                                 "The alignment can't have no matches here");

            ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleDeletion done: " << headAlignment << "-" << tailAlignment);
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleDeletion having gap is worse than keeping the originals: " <<
                                        headAlignment.smithWatermanScore << "<=" << sws << ":" << headMismatches << "+" << bestMismatches);
        }
    }
    else
    {
        ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleDeletion gap has too many mismatches at the flanks: bestLeftFlankMismatches=" <<
                                    bestLeftFlankMismatches << " bestRightFlankMismatches=" << bestRightFlankMismatches << " bestOffset=" << bestOffset);
    }
}


/**
 * \brief Patches the front fragment with cigar that produces the lowest number of mismatches assuming there
 *        is an insertion in the read somewhere between the headAlignment first seed and tailAlignment first seed
 *
 * \headAlignment alignment where the leftmost part matches the reference
 * \param headSeedOffset offset of the seed for headAlignment from the left-most base in respect to the reference
 * \tailAlignment alignment where the rightmost part matches the reference
 * \param tailSeedOffset offset of the seed for tailAlignment from the left-most base in respect to the reference
 */
void SimpleIndelAligner::alignSimpleInsertion(
    Cigar &cigarBuffer,
    FragmentMetadata &headAlignment,
    const unsigned headSeedOffset,
    const unsigned headSeedLength,
    FragmentMetadata &tailAlignment,
    const unsigned tailSeedOffset,
    const unsigned tailSeedLength,
    const std::vector<reference::Contig> &contigList,
    const flowcell::ReadMetadataList &readMetadataList) const
{
    if (headSeedOffset < headAlignment.getBeginClippedLength())
    {
        ISAAC_THREAD_CERR_DEV_TRACE("alignSimpleInsertion: clipping on the left overlaps the head seed.");
        return;
    }

    if (tailAlignment.getBeginClippedLength() + tailAlignment.getObservedLength() < tailSeedOffset + tailSeedLength)
    {
        ISAAC_THREAD_CERR_DEV_TRACE("alignSimpleInsertion: clipping on the right overlaps the tail seed.");
        return;
    }

    ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleInsertion: " << headAlignment << "-" << tailAlignment);

    // Don't allow insertion to be placed within the first seed to capture the earliest possible location as
    // insertions directly reduce number of mismatches, thus causing unfair competition against non-gapped alignment candidates
    const unsigned tailOffset = headSeedOffset + headSeedLength;

    const unsigned observedEnd = tailAlignment.getBeginClippedLength() + tailAlignment.getObservedLength();
    const unsigned insertionLength = boost::numeric_cast<unsigned>(headAlignment.getUnclippedPosition() - tailAlignment.getUnclippedPosition());

    if (tailSeedOffset - headSeedOffset < insertionLength + headSeedLength)
    {
        ISAAC_THREAD_CERR_DEV_TRACE("alignSimpleInsertion: insertion too long (must fit between the anchoring seeds):" << insertionLength << " max: " << (tailSeedOffset - headSeedOffset));
        return;
    }

    const Read &read = headAlignment.getRead();
    const bool reverse = headAlignment.reverse;
    const std::vector<char>::const_iterator sequenceBegin = read.getStrandSequence(reverse).begin();
    const std::vector<char> &reference = contigList[headAlignment.contigId].forward_;

    std::vector<char>::const_iterator tailIterator = sequenceBegin + tailOffset + insertionLength;
    unsigned tailLength = observedEnd - tailOffset - insertionLength;
    // number of mismatches when insertion is at the left extremity
    ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleInsertion insertionLength:" << insertionLength << " tailLength:" << tailLength);
    const unsigned tailMismatches = countMismatches(tailIterator,
                                                    reference.begin() + headAlignment.getUnclippedPosition() + tailOffset, reference.end(),
                                                    tailLength, &boost::cref<char>);

    unsigned leftFlankMismatches = countMismatches(tailIterator - insertionLength - GAP_FLANK_BASES,
                                                   reference.begin() + headAlignment.getUnclippedPosition() + tailOffset - GAP_FLANK_BASES, reference.end(),
                                                   GAP_FLANK_BASES, &boost::cref<char>);;

    unsigned rightFlankMismatches = countMismatches(tailIterator,
                                                    reference.begin() + headAlignment.getUnclippedPosition() + tailOffset, reference.end(),
                                                    std::min(GAP_FLANK_BASES, tailLength), &boost::cref<char>);;

    ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleInsertion " <<
                                "sequence:" << std::string(read.getStrandSequence(reverse).begin(), read.getStrandSequence(reverse).end()) <<
                                " ref1=" << std::string(reference.begin() + headAlignment.position,
                                                        std::min(reference.end(), reference.begin() + headAlignment.position + headAlignment.getObservedLength())) <<
                                " ref2=" << std::string(reference.begin() + tailAlignment.position,
                                                        std::min(reference.end(), reference.begin() + tailAlignment.position + tailAlignment.getObservedLength())));
/*
    if (!tailMismatches)
    {
        ISAAC_THREAD_CERR_DEV_TRACE("alignSimpleInsertion: no point to try, the head alignment is already good enough");
        return;
    }*/


    // try to introduce the insertion see if the number of mismatches reduces below the original

    // number of mismatches when insertion is at the leftmost possible position
    unsigned rightRealignedMismatches = tailMismatches;/*countMismatches(tailIterator,
                                                        reference.begin() + headAlignment.position + tailOffset, reference.end(),
                                                        tailLength, &boost::cref<char>);*/
    // we're starting at the situation where the whole tail of the head alignment is moved by -insertionLength
    unsigned leftRealignedMismatches = 0;

    ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleInsertion" <<
                                " insertionLength=" << insertionLength <<
                                " tailLength=" << tailLength << "tail:" << std::string(tailIterator, tailIterator + tailLength) <<
                                " tailMismatches=" << tailMismatches <<
                                " rightRealignedMismatches=" << rightRealignedMismatches <<
                                " leftRealignedMismatches=" << leftRealignedMismatches);

    std::vector<char>::const_iterator referenceIterator = reference.begin() + headAlignment.getUnclippedPosition() + tailOffset;

    unsigned bestMismatches = tailMismatches;
    unsigned bestOffset = tailOffset;
    unsigned bestLeftFlankMismatches = leftFlankMismatches;
    unsigned bestRightFlankMismatches = rightFlankMismatches;

    ISAAC_THREAD_CERR_DEV_TRACE("headTailOffset=" << tailOffset << " tailSeedOffset=" << tailSeedOffset << " insertionLength=" << insertionLength);

    for (unsigned insertionOffset = tailOffset; bestMismatches && insertionOffset <= tailSeedOffset - insertionLength;
        ++insertionOffset, ++tailIterator, ++referenceIterator, --tailLength)
    {
        const unsigned thisOffsetMismatches = leftRealignedMismatches + rightRealignedMismatches;

        if (bestMismatches > thisOffsetMismatches)
        {
            bestOffset = insertionOffset;
            bestMismatches = thisOffsetMismatches;
            bestLeftFlankMismatches = leftFlankMismatches;
            bestRightFlankMismatches = rightFlankMismatches;
        }
        ISAAC_THREAD_CERR_DEV_TRACE("bestMismatches=" << bestMismatches <<
                                    " leftRealignedMismatches=" << leftRealignedMismatches <<
                                    " rightRealignedMismatches=" << rightRealignedMismatches << " insertionOffset=" << insertionOffset <<
                                    " leftFlankMismatches" << leftFlankMismatches <<
                                    " rightFlankMismatches" << rightFlankMismatches);
//        ISAAC_THREAD_CERR_DEV_TRACE(" headTail=" << std::string(tailIterator, read.getStrandSequence(reverse).end()) <<
//                                    " headTailRef=" << std::string(referenceIterator, referenceIterator + std::distance(tailIterator, read.getStrandSequence(reverse).end())));

        // unlike deletions, insertions consume some bases of the sequence. Thus, we check the base that exits the insertion on the left side,
        const bool newLeftMismatch = !isMatch(*(tailIterator - insertionLength), *referenceIterator);
        leftRealignedMismatches += newLeftMismatch;
        leftFlankMismatches += newLeftMismatch;
        if (insertionOffset >= GAP_FLANK_BASES)
        {
            leftFlankMismatches -= !isMatch(*(tailIterator - insertionLength - GAP_FLANK_BASES), *(referenceIterator - GAP_FLANK_BASES));
        }
        // and the one that enters it on the right side
        const bool disappearingRightMismatch = !isMatch(*tailIterator, *(referenceIterator));
        rightRealignedMismatches -= disappearingRightMismatch;
        rightFlankMismatches -= disappearingRightMismatch;
        if (tailLength > GAP_FLANK_BASES)
        {
            rightFlankMismatches += !isMatch(*(tailIterator + GAP_FLANK_BASES), *(referenceIterator + GAP_FLANK_BASES));
        }


    }

    const long clippingPositionOffset = headAlignment.getBeginClippedLength();
    const unsigned leftMapped = bestOffset - clippingPositionOffset;
    ISAAC_ASSERT_MSG(leftMapped, "Simple insertions are not allowed to be placed at the very beginning of the read")
    const unsigned headMismatches = countMismatches(sequenceBegin + clippingPositionOffset,
                                                    reference.begin() + headAlignment.position, reference.end(),
                                                    leftMapped, &boost::cref<char>);

    const unsigned newMismatches = headMismatches + bestMismatches;
    const unsigned sws = normalizedMismatchScore_ * newMismatches + normalizedGapOpenScore_ +
        std::min(normalizedMaxGapExtendScore_, (insertionLength - 1) * normalizedGapExtendScore_);
    if (bestLeftFlankMismatches <= GAP_FLANK_MISMATCHES_MAX && bestRightFlankMismatches <= GAP_FLANK_MISMATCHES_MAX)
    {
        if (tailAlignment.smithWatermanScore > sws || (tailAlignment.smithWatermanScore == sws && tailAlignment.getMismatchCount() > newMismatches))
        {
            const unsigned cigarOffset = cigarBuffer.size();


            if (clippingPositionOffset)
            {
                cigarBuffer.addOperation(clippingPositionOffset, Cigar::SOFT_CLIP);
            }

            cigarBuffer.addOperation(leftMapped, Cigar::ALIGN);
            cigarBuffer.addOperation(insertionLength, Cigar::INSERT);

            const unsigned rightMapped = headAlignment.getObservedLength()  + headAlignment.getEndClippedLength() - leftMapped - tailAlignment.getEndClippedLength() - insertionLength;
            ISAAC_ASSERT_MSG(rightMapped, "Simple insertions are not allowed to be placed at the very end of the read " <<
                             headAlignment << "-" << tailAlignment << " leftMapped:" << leftMapped << " insertionLength:" << insertionLength <<
                             " tailSeedOffset:" << tailSeedOffset << " headSeedOffset:" << headSeedOffset);

            cigarBuffer.addOperation(rightMapped, Cigar::ALIGN);

            const unsigned clipEndBases = tailAlignment.getEndClippedLength();
            if (clipEndBases)
            {
                cigarBuffer.addOperation(clipEndBases, Cigar::SOFT_CLIP);
            }

            // Unlike alignSimpleDeletion, update tailAlignment as it is ordered earlier than headAlignment in the fragment list
            tailAlignment.resetAlignment(cigarBuffer);
            // carry over alignment-independent clipping information (quality trimming, adapter masking)
            tailAlignment.leftClipped() = headAlignment.leftClipped();
            ISAAC_ASSERT_MSG(updateFragmentCigar(readMetadataList, reference, tailAlignment,
                                                 headAlignment.position, cigarBuffer, cigarOffset),
                                                 "The alignment can't have no matches here");

            ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleInsertion done: " << headAlignment << "-" << tailAlignment);
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleInsertion having gap is worse than keeping the originals: " <<
                                        tailAlignment.smithWatermanScore << "<=" << sws << ":" << headMismatches << "+" << bestMismatches);
        }
    }
    else
    {
        ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleInsertion gap has too many mismatches at the flanks: bestLeftFlankMismatches=" <<
                                    bestLeftFlankMismatches << " bestRightFlankMismatches=" << bestRightFlankMismatches << " bestOffset=" << bestOffset);
    }
}

/**
 ** \brief Comparison of FragmentMetadata by unclipped reference position
 **/
bool orderByUnclippedPosition(const FragmentMetadata &left, const FragmentMetadata &right)
{
    return
        (left.contigId < right.contigId || (left.contigId == right.contigId &&
            (left.getUnclippedPosition() < right.getUnclippedPosition()))
        );
}


/**
 * \brief Catches the cases of single indel in the fragment by analyzing the
 *        seed alignment conflicts.
 *
 * \precondition The FragmentMetadata::firstSeedIndex is valid
 * \precondition The fragmentList is ordered
 * \precondition The fragmentList contains unique alignments only
 */
void SimpleIndelAligner::alignSimpleIndels(
    Cigar &cigarBuffer,
    const std::vector<reference::Contig> &contigList,
    const flowcell::ReadMetadataList &readMetadataList,
    const SeedMetadataList &seedMetadataList,
    FragmentMetadataList &fragmentList) const
{
    if (fragmentList.size() < 2)
    {
        return;
    }

    std::sort(fragmentList.begin(), fragmentList.end(), orderByUnclippedPosition);

    FragmentMetadataList::iterator head = fragmentList.begin();
    // evaluate pairs of alignments (they are ordered). If the indel alignment is better than any of the members of the pair,
    // patch the front (lowest alignment position) one.
    // Use the back one for evaluation of the subsequent pair.
    for (std::vector<FragmentMetadata>::iterator tail = head + 1; fragmentList.end() != tail; ++tail, ++head)
    {
        if (head->contigId == tail->contigId && head->reverse == tail->reverse)
        {
            const SeedMetadata &headSeed = seedMetadataList.at(head->firstSeedIndex);
            const SeedMetadata &tailSeed = seedMetadataList.at(tail->firstSeedIndex);

            // ignore conflicting repeat alignments.
            //if(headSeed.getIndex() != tailSeed.getIndex())
            {
                const long distance = tail->getUnclippedPosition() - head->getUnclippedPosition();
                ISAAC_ASSERT_MSG(distance, "distance must be non-zero for gap introduction");
                if (std::abs(distance) < semialignedGapLimit_)
                {
                    ISAAC_ASSERT_MSG(0 <= distance, "Pre-condition violation. The alignments are expected to be unique and ordered. Got: " << *head << "-" << *tail);
                    const long headSeedOffset = head->reverse ? head->getReadLength() - headSeed.getOffset() - headSeed.getLength() : headSeed.getOffset();
                    const long tailSeedOffset = head->reverse ? head->getReadLength() - tailSeed.getOffset() - tailSeed.getLength() : tailSeed.getOffset();
                    const long expectedSeedDistance = tailSeedOffset - headSeedOffset;
                    ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleIndels " << headSeed << headSeedOffset << "headSeed, " << tailSeed << tailSeedOffset << "tailSeed, expectedSeedDistance:" << expectedSeedDistance);
                    if (0 < expectedSeedDistance)
                    {
                        // seeds are ordered same way as the alignments. First alignment is before the second one.
                        // This is a deletion in the data with distance being the length of the deletion.
                        alignSimpleDeletion(cigarBuffer, *head, headSeedOffset, *tail,
                                            tailSeedOffset, tailSeed.getLength(), contigList, readMetadataList);
                    }
                    else
                    {
                        // seeds are ordered in the opposite to alignments. This means they are closer than they should be.
                        alignSimpleInsertion(cigarBuffer, *tail, tailSeedOffset, tailSeed.getLength(), *head,
                                             headSeedOffset, headSeed.getLength(), contigList, readMetadataList);
                    }
                }
                else
                {
                    ISAAC_THREAD_CERR_DEV_TRACE(" alignSimpleIndels " << headSeed << "headSeed " << tailSeed << "tailSeed too far apart: " << distance);
                }
            }
        }
    }
}


} // namespace fragmentBuilder
} // namespace alignment
} // namespace isaac
