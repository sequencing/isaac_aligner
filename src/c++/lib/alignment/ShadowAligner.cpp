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
 ** \file ShadowAligner.cpp
 **
 ** \brief See ShadowAligned.hh
 ** 
 ** \author Come Raczy
 **/

#include <algorithm>
#include <boost/format.hpp>

#include "alignment/Quality.hh"
#include "alignment/ShadowAligner.hh"
#include "alignment/matchSelector/FragmentSequencingAdapterClipper.hh"
#include "oligo/Kmer.hh"
#include "oligo/KmerGenerator.hpp"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

ShadowAligner::ShadowAligner(const flowcell::FlowcellLayoutList &flowcellLayoutList,
                             const unsigned gappedMismatchesMax,
                             const bool avoidSmithWaterman,
                             const int gapMatchScore,
                             const int gapMismatchScore,
                             const int gapOpenScore,
                             const int gapExtendScore,
                             const int minGapExtendScore)
    : gappedMismatchesMax_(gappedMismatchesMax),
      ungappedAligner_(gapMatchScore, gapMismatchScore, gapOpenScore, gapExtendScore, minGapExtendScore),
      gappedAligner_(flowcellLayoutList, avoidSmithWaterman, gapMatchScore, gapMismatchScore, gapOpenScore, gapExtendScore, minGapExtendScore),
      shadowCigarBuffer_(Cigar::getMaxOperationsForReads(flowcellLayoutList) *
                         unreasonablyHighDifferenceBetweenMaxAndMinInsertSizePlusFlanks_)
{
    shadowCandidatePositions_.reserve(unreasonablyHighDifferenceBetweenMaxAndMinInsertSizePlusFlanks_);
    shadowKmerPositions_.reserve(shadowKmerCount_);
}

unsigned ShadowAligner::hashShadowKmers(const std::vector<char> &sequence)
{
    shadowKmerPositions_.clear();
    // initialize all k-mers to the magic value -1 (NOT_FOUND)
    shadowKmerPositions_.resize(shadowKmerCount_, -1);
    // 
    oligo::KmerGenerator<unsigned> kmerGenerator(sequence.begin(), sequence.end(), shadowKmerLength_);
    unsigned positionsCount = 0;
    unsigned kmer;
    std::vector<char>::const_iterator position;
    while (kmerGenerator.next(kmer, position))
    {
        if (-1 == shadowKmerPositions_[kmer])
        {
            shadowKmerPositions_[kmer] = (position - sequence.begin());
            ++positionsCount;
        }
    }
    return positionsCount;
}

void ShadowAligner::findShadowCandidatePositions(
    const std::vector<char>::const_iterator referenceBegin,
    const std::vector<char>::const_iterator referenceEnd,
    const std::vector<char> &shadowSequence)
{
    hashShadowKmers(shadowSequence);

    // find matching positions in the reference by k-mer comparison
    oligo::KmerGenerator<unsigned> kmerGenerator(referenceBegin, referenceEnd, shadowKmerLength_);
    unsigned kmer;
    std::vector<char>::const_iterator position;
    while(kmerGenerator.next(kmer, position))
    {
        if (-1 != shadowKmerPositions_[kmer])
        {
            const long candidatePosition = position - referenceBegin - shadowKmerPositions_[kmer];
            // avoid spurious repetitions of start positions
            if (shadowCandidatePositions_.empty() || shadowCandidatePositions_.back() != candidatePosition)
            {
                if (shadowCandidatePositions_.size() == shadowCandidatePositions_.capacity())
                {
                    // too many candidate positions. Just stop here. The alignment score will be miserable anyway.
                    break;
                }
                shadowCandidatePositions_.push_back(candidatePosition);
//                ISAAC_THREAD_CERR << "kmer: " << oligo::Bases<2, oligo::Kmer>(kmer, shadowKmerLength_) << " " << candidatePosition << std::endl;
            }
        }
    }

    // remove duplicate positions
    if (!shadowCandidatePositions_.empty())
    {
        std::sort(shadowCandidatePositions_.begin(), shadowCandidatePositions_.end());
        shadowCandidatePositions_.erase(std::unique(shadowCandidatePositions_.begin(),
                                                    shadowCandidatePositions_.end()),
                                        shadowCandidatePositions_.end());
    }
}

/**
 * \brief if the best template is longer than the dominant template, attempt to rescue shadow
 *        within a wider range
 * \param bestTemplateLength 0 indicates there is no best template
 */
std::pair <long, long> calculateShadowRescueRange(
    const FragmentMetadata &orphan,
    const TemplateLengthStatistics &templateLengthStatistics,
    const long bestTemplateLength)
{
    const Cluster &cluster = orphan.getCluster();

    const unsigned shadowReadIndex = (orphan.readIndex + 1) % 2;
    const unsigned readLengths[] = {cluster[0].getLength(), cluster[1].getLength()};
    long shadowMinPosition = templateLengthStatistics.mateMinPosition(orphan.readIndex, orphan.reverse, orphan.position, readLengths);
    long shadowMaxPosition = templateLengthStatistics.mateMaxPosition(orphan.readIndex, orphan.reverse, orphan.position, readLengths) +
        readLengths[shadowReadIndex] - 1;

    if (bestTemplateLength)
    {
        if (shadowMinPosition < long(orphan.getFStrandReferencePosition().getPosition()))
        {
            shadowMinPosition = std::min(long(orphan.getRStrandReferencePosition().getPosition()) - bestTemplateLength, shadowMinPosition);
        }

        // Notice we're not using orphan.getRStrandReferencePosition in the comparison as orphans with long deletions have
        // their end often well outside the dominant template length
        if (shadowMaxPosition > long(orphan.getFStrandReferencePosition().getPosition()))
        {
            shadowMaxPosition = std::max(long(orphan.getFStrandReferencePosition().getPosition()) + bestTemplateLength, shadowMaxPosition);
        }
    }

    const std::pair <long, long> ret(shadowMinPosition - 10, shadowMaxPosition + 10);
    return ret;
}

/**
 * \return false when no reasonable placement for shadow found. If at that point the shadowList is not empty,
 * this means that the shadow falls at a repetitive region and rescuing should not be considered
 */
bool ShadowAligner::rescueShadow(
    const std::vector<reference::Contig> &contigList,
    const FragmentMetadata &orphan,
    std::vector<FragmentMetadata> &shadowList,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const TemplateLengthStatistics &templateLengthStatistics,
    const long bestTemplateLength)
{
    if (!templateLengthStatistics.isCoherent())
    {
        ISAAC_THREAD_CERR_DEV_TRACE("    Rescuing impossible. Incoherent tls");
        return false;
    }
    shadowCigarBuffer_.clear();
    ISAAC_ASSERT_MSG(2 > orphan.readIndex, "Paired reads means 2");
    const Cluster &cluster = orphan.getCluster();
    const unsigned shadowReadIndex = (orphan.readIndex + 1) % 2;
    const Read &shadowRead = cluster[shadowReadIndex];
    // identify the orientation and range of reference positions of the orphan
    const reference::Contig &contig = contigList[orphan.contigId];
    const bool shadowReverse = templateLengthStatistics.mateOrientation(orphan.readIndex, orphan.reverse);
    const std::vector<char> &reference = contig.forward_;
    const std::pair<long, long> shadowRescueRange = calculateShadowRescueRange(orphan, templateLengthStatistics, bestTemplateLength);
    if (shadowRescueRange.second < shadowRescueRange.first)
    {
        ISAAC_THREAD_CERR_DEV_TRACE((boost::format("    Rescuing impossible: shadowMaxPosition < shadowMinPosition "
            "%l < %l") % shadowRescueRange.second % shadowRescueRange.first).str());
        return false;
    }
    if (shadowRescueRange.second + 1 + shadowRead.getLength() < 0)
    {
        ISAAC_THREAD_CERR_DEV_TRACE((boost::format("    Rescuing impossible: shadowMaxPosition + 1 + shadowRead.getLength() < 0 "
            "%l < 0") % (shadowRescueRange.second + 1 + shadowRead.getLength())).str());
        return false;
    }
    // find all the candidate positions for the shadow on the identified reference region
    shadowCandidatePositions_.clear();
    const std::vector<char> &shadowSequence = shadowReverse ? shadowRead.getReverseSequence() : shadowRead.getForwardSequence();
    const long candidatePositionOffset = std::max(0L, shadowRescueRange.first);
    findShadowCandidatePositions(
        reference.begin() + candidatePositionOffset,
        reference.begin() + std::min((long)reference.size(), shadowRescueRange.second + 1),
        shadowSequence);

    ISAAC_THREAD_CERR_DEV_TRACE("findShadowCandidatePositions found " << shadowCandidatePositions_.size() << " positions in range [" <<
                                (candidatePositionOffset) << ";" <<
                                (std::min((long)reference.size(), shadowRescueRange.second + 1)) << "]" <<
                                " bestTemplateLength:" << bestTemplateLength);
    // align the shadow to the candidate positions and keep the best fragment
    shadowList.clear();

    matchSelector::FragmentSequencingAdapterClipper adapterClipper(sequencingAdapters);

    FragmentMetadata *bestFragment = 0;
    BOOST_FOREACH(long strandPosition, shadowCandidatePositions_)
    {
        if (shadowList.size() == shadowList.capacity())
        {
            return false;
        }
        strandPosition += candidatePositionOffset;
        FragmentMetadata fragment(&cluster, &shadowCigarBuffer_, shadowReadIndex);
        fragment.reverse = shadowReverse;
        fragment.contigId = orphan.contigId;
        fragment.position = strandPosition;

        adapterClipper.checkInitStrand(fragment, contig);
        if (ungappedAligner_.alignUngapped(fragment, shadowCigarBuffer_, readMetadataList, adapterClipper, contig))
        {
            ISAAC_THREAD_CERR_DEV_TRACE("    Rescuing: Aligned: " << fragment);
            shadowList.push_back(fragment);
            if (0 == bestFragment || ISAAC_LP_LESS(bestFragment->logProbability, fragment.logProbability))
            {
                bestFragment = &shadowList.back();
            }
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE("    Rescuing: Unaligned: " << fragment);
        }
    }

    if (!bestFragment)
    {
        return false;
    }

    if(BandedSmithWaterman::mismatchesCutoff < bestFragment->mismatchCount)
    {
        FragmentMetadataList::const_iterator nextCandidate = shadowList.begin();
        BOOST_FOREACH(FragmentMetadata &fragment, shadowList)
        {
            ++nextCandidate;
            if (shadowList.end() != nextCandidate &&
                nextCandidate->position - fragment.position < BandedSmithWaterman::distanceCutoff)
            {
                // Use the gapped aligner if necessary
                if (BandedSmithWaterman::mismatchesCutoff < fragment.mismatchCount)
                {
                    FragmentMetadata tmp = fragment;
                    const unsigned matchCount = gappedAligner_.alignGapped(tmp, shadowCigarBuffer_,
                                                                             readMetadataList, adapterClipper, contig);
                    ISAAC_THREAD_CERR_DEV_TRACE("    Rescuing:     Gap-aligned: " << tmp);
                    if (matchCount && matchCount + BandedSmithWaterman::WIDEST_GAP_SIZE > fragment.getObservedLength() &&
                        (tmp.mismatchCount <= gappedMismatchesMax_) &&
                        (fragment.mismatchCount > tmp.mismatchCount) &&
                        ISAAC_LP_LESS(fragment.logProbability, tmp.logProbability))
                    {
                        fragment = tmp;
                        if (ISAAC_LP_LESS(bestFragment->logProbability, fragment.logProbability))
                        {
                            bestFragment = &fragment;
                        }
                    }
                }
            }
        }
    }

/*    if(bestFragment->mismatchCount > bestFragment->getReadLength() / 8 ||
        bestFragment->logProbability < LOG_MISMATCH_Q40 / 4 * bestFragment->getReadLength())
    {
        ISAAC_THREAD_CERR_DEV_TRACE("    Rescued shadow too bad: " << *bestFragment);

        shadowList.clear();
        return false;
    }
    else*/
    {
        if (&shadowList.front() != bestFragment)
        {
            std::swap(shadowList.front(), *bestFragment);
        }
        return true;
    }
}

} // namespace alignment
} // namespace isaac
