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
#include "oligo/KmerGenerator.hpp"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

ShadowAligner::ShadowAligner(const flowcell::FlowcellLayoutList &flowcellLayoutList,
                             const unsigned gappedMismatchesMax,
                             const FragmentBuilder &fragmentBuilder)
    : gappedMismatchesMax_(gappedMismatchesMax),
      fragmentBuilder_(fragmentBuilder),
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
    const unsigned int kmers = hashShadowKmers(shadowSequence);
    if (20 > kmers)
    {
        ISAAC_THREAD_CERR_DEV_TRACE("findShadowCandidatePositions shadow is not diverse enough. Found " << kmers << " " << shadowKmerLength_ << "-mers");
        // shadow is not diverse enough
        return;
    }
    // find matching positions in the reference by k-mer comparison
    oligo::KmerGenerator<unsigned> kmerGenerator(referenceBegin, referenceEnd, shadowKmerLength_);
    unsigned kmer;
    std::vector<char>::const_iterator position;
    while(kmerGenerator.next(kmer, position))
    {
        if (-1 != shadowKmerPositions_[kmer])
        {
            const long candidatePosition = position - referenceBegin - shadowKmerPositions_[kmer];
            // avoid spurions repetitions of start positions
            if (shadowCandidatePositions_.empty() || shadowCandidatePositions_.back() != candidatePosition)
            {
                shadowCandidatePositions_.push_back(candidatePosition);
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
 * \return false when no reasonable placement for shadow found. If at that point the shadowList is not empty,
 * this means that the shadow falls at a repetitive region and rescuing should not be considered
 */
bool ShadowAligner::rescueShadow(
    const std::vector<reference::Contig> &contigList,
    const FragmentMetadata &orphan,
    std::vector<FragmentMetadata> &shadowList,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const TemplateLengthStatistics &templateLengthStatistics)
{
    if (!templateLengthStatistics.isCoherent())
    {
        ISAAC_THREAD_CERR_DEV_TRACE("    Rescuing impossible. Incoherent tls");
        return false;
    }
    shadowCigarBuffer_.clear();
    assert(2 > orphan.readIndex);
    const Cluster &cluster = orphan.getCluster();
    const unsigned shadowReadIndex = (orphan.readIndex + 1) % 2;
    const Read &shadowRead = cluster[shadowReadIndex];
    // identify the orientation and range of reference positions of the orphan
    const unsigned readLengths[] = {cluster[0].getLength(), cluster[1].getLength()};
    const reference::Contig &contig = contigList[orphan.contigId];
    const bool shadowReverse = templateLengthStatistics.mateOrientation(orphan.readIndex, orphan.reverse);
    const std::vector<char> &reference = contig.forward_;
    const long shadowMinPosition = templateLengthStatistics.mateMinPosition(orphan.readIndex, orphan.reverse, orphan.position, readLengths) - 10;
    const long shadowMaxPosition = templateLengthStatistics.mateMaxPosition(orphan.readIndex, orphan.reverse, orphan.position, readLengths) + 
        cluster[shadowReadIndex].getLength() - 1 + 10;
    if (shadowMaxPosition < shadowMinPosition)
    {
        ISAAC_THREAD_CERR_DEV_TRACE((boost::format("    Rescuing impossible: shadowMaxPosition < shadowMinPosition "
            "%l < %l") % shadowMaxPosition % shadowMinPosition).str());
        return false;
    }
    if (shadowMaxPosition + 1 + shadowRead.getLength() < 0)
    {
        ISAAC_THREAD_CERR_DEV_TRACE((boost::format("    Rescuing impossible: shadowMaxPosition + 1 + shadowRead.getLength() < 0 "
            "%l < 0") % (shadowMaxPosition + 1 + shadowRead.getLength())).str());
        return false;
    }
    // find all the candidate positions for the shadow on the identified reference region
    shadowCandidatePositions_.clear();
    const std::vector<char> &shadowSequence = shadowReverse ? shadowRead.getReverseSequence() : shadowRead.getForwardSequence();
    const long candidatePositionOffset = std::max(0L, shadowMinPosition);
    findShadowCandidatePositions(
        reference.begin() + candidatePositionOffset,
        reference.begin() + std::min((long)reference.size(), shadowMaxPosition + 1),
        shadowSequence);

    ISAAC_THREAD_CERR_DEV_TRACE("findShadowCandidatePositions found " << shadowCandidatePositions_.size() << " positions in range [" <<
                                (candidatePositionOffset) << ";" <<
                                (std::min((long)reference.size(), shadowMaxPosition + 1)) << "]");
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
        fragmentBuilder_.alignUngapped(fragment, shadowCigarBuffer_, readMetadataList, adapterClipper, contig);
        using boost::format;
        ISAAC_THREAD_CERR_DEV_TRACE("    Rescuing: Aligned: " << fragment);
        shadowList.push_back(fragment);
        if (0 == bestFragment || ISAAC_LP_LESS(bestFragment->logProbability, fragment.logProbability))
        {
            bestFragment = &shadowList.back();
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
                    const unsigned matchCount = fragmentBuilder_.alignGapped(tmp, shadowCigarBuffer_,
                                                                             readMetadataList, adapterClipper, contig);
                    ISAAC_THREAD_CERR_DEV_TRACE((boost::format("    Rescuing:     Gap-aligned: %s") % tmp).str());
                    if (matchCount && matchCount + BandedSmithWaterman::widestGapSize > fragment.getObservedLength() &&
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
