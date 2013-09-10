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
 ** \file FragmentMetadata.hh
 **
 ** \brief Component to encapsulate all metadata associated to a fragment.
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_METADATA_HH
#define iSAAC_ALIGNMENT_FRAGMENT_METADATA_HH

#include <vector>
#include <algorithm>
#include <iostream>

#include "alignment/Alignment.hh"
#include "alignment/Cigar.hh"
#include "alignment/Cluster.hh"
#include "alignment/SeedId.hh"
#include "common/FiniteCapacityVector.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Alignment information for a fragment (as defined by the SAM Format
 ** Specification (v1.4-r962)
 **
 ** This component is the building block for the FragmentBuilder. It is
 ** designed for efficiency, and does not involve any memory allocation (which
 ** is the reason why it does not and should not store the CIGAR).
 **/
struct FragmentMetadata
{
    friend std::ostream &operator<<(std::ostream &os, const FragmentMetadata &f);

    /**
     * \param readLength is needed to preallocate buffer to avoid memory operations during the processing.
     *        auxiliary code, such as unit tests, does not need to supply it.
     */
    FragmentMetadata():
        cluster(0), contigId(reference::ReferencePosition::MAX_CONTIG_ID), position(0), lowClipped(0), highClipped(0),
        observedLength(0), readIndex(0), reverse(false), cigarOffset(0),
        cigarLength(0), cigarBuffer(0), mismatchCount(0), matchesInARow(0), gapCount(0), editDistance(0), logProbability(0.0),
        firstSeedIndex(-1), repeatSeedsCount(0), uniqueSeedCount(0), nonUniqueSeedOffsets(std::numeric_limits<unsigned>::max(), 0UL),
        alignmentScore(-1U),
        smithWatermanScore(0)
    {
        std::fill(mismatchCycles, mismatchCycles + maxCycles_, 0);
    }

    FragmentMetadata(const Cluster *cluster, const std::vector<unsigned> *cigarBuffer, unsigned readIndex):
        cluster(cluster), contigId(reference::ReferencePosition::MAX_CONTIG_ID), position(0), lowClipped(0), highClipped(0),
        observedLength(0), readIndex(readIndex), reverse(false), cigarOffset(0),
        cigarLength(0), cigarBuffer(cigarBuffer), mismatchCount(0), matchesInARow(0), gapCount(0), editDistance(0), logProbability(0.0),
        firstSeedIndex(-1), repeatSeedsCount(0), uniqueSeedCount(0), nonUniqueSeedOffsets(std::numeric_limits<unsigned>::max(), 0UL),
        alignmentScore(-1U),
        smithWatermanScore(0)
    {
        std::fill(mismatchCycles, mismatchCycles + maxCycles_, 0);
    }

    bool isReverse() const {return reverse;}
    unsigned getReadLength() const
    {
        assert(0 != cluster);
        assert(cluster->size() > readIndex);
        return (*cluster)[readIndex].getLength();
    }
    unsigned getReadIndex() const {return readIndex;}
    unsigned getObservedLength() const {return isAligned() ? observedLength : 0;}
    unsigned getAlignmentScore() const {return alignmentScore;}
    void setAlignmentScore(unsigned as) {alignmentScore = as;}
    unsigned getCigarLength() const {return cigarLength;}
    /// Position of the first base of the fragment
    reference::ReferencePosition getFStrandReferencePosition() const
    {
        return !isNoMatch() ?
            reference::ReferencePosition(contigId, position) :
            reference::ReferencePosition(reference::ReferencePosition::NoMatch);
    }
    /// Position of the last base of the fragment
    reference::ReferencePosition getRStrandReferencePosition() const
    {
        // ensure that the position is positive!
        return !isNoMatch() ?
            reference::ReferencePosition(contigId, std::max(position + observedLength, 1L) - 1) :
            reference::ReferencePosition(reference::ReferencePosition::NoMatch);
    }
    /// Position of the fragment
    reference::ReferencePosition getStrandReferencePosition() const {
        return isReverse() ? getRStrandReferencePosition() : getFStrandReferencePosition();
    }

    /// Same as f-strand position
    reference::ReferencePosition getBeginReferencePosition() const
    {
        return getFStrandReferencePosition();
    }

    /// Different from r-strand position in that it always points to the base following the last unclipped base of the fragment
    reference::ReferencePosition getEndReferencePosition() const
    {
        return !isNoMatch() ?
            reference::ReferencePosition(contigId, position + observedLength) :
            reference::ReferencePosition(reference::ReferencePosition::NoMatch);
    }

/// First cycle of fragment bcl data
    std::vector<char>::const_iterator getBclData() const {
        return cluster->getBclData(getReadIndex());
    }
    /// Cluster of the fragment
    const Cluster &getCluster() const {
        return *cluster;
    }

    const Read &getRead() const {
        return getCluster()[getReadIndex()];
    }

    Cigar::const_iterator cigarBegin() const
    {
        ISAAC_ASSERT_MSG(isAligned(), "Requesting CIGAR of unaligned fragment is not allowed");
        return cigarBuffer->begin() + cigarOffset;
    }

    Cigar::const_iterator cigarEnd() const
    {
        ISAAC_ASSERT_MSG(isAligned(), "Requesting CIGAR of unaligned fragment is not allowed");
        return cigarBuffer->begin() + cigarOffset + cigarLength;
    }

    long getBeginClippedLength() const
    {
        if (cigarBuffer && cigarLength)
        {
            std::pair<unsigned, Cigar::OpCode> operation = Cigar::decode(cigarBuffer->at(cigarOffset));
            if (Cigar::SOFT_CLIP == operation.second)
            {
                return operation.first;
            }
        }
        return 0;
    }

    long getEndClippedLength() const
    {
        if (cigarBuffer && cigarLength)
        {
            std::pair<unsigned, Cigar::OpCode> operation = Cigar::decode(cigarBuffer->at(cigarOffset + cigarLength - 1));
            if (Cigar::SOFT_CLIP == operation.second)
            {
                return operation.first;
            }
        }
        return 0;
    }

    /// Unlike the observed length, excludes gaps (deletions and insertion bases)
    unsigned getMappedLength() const
    {
        ISAAC_ASSERT_MSG(cigarBuffer && cigarLength, "Read must have a valid CIGAR");
        return Cigar::getMappedLength(cigarBuffer->begin() + cigarOffset,
                                      cigarBuffer->begin() + cigarOffset + cigarLength);
    }

    /**
     * \brief Returns unadjusted position if it is adjusted due to a soft clipping
     */
    long getUnclippedPosition() const
    {
        return position - getBeginClippedLength();
    }

    unsigned getMismatchCount() const {
        return mismatchCount;
    }

    unsigned getGapCount() const {
        return gapCount;
    }

    unsigned getEditDistance() const {
        return editDistance;
    }

    const unsigned short* getMismatchCyclesBegin() const {return mismatchCycles;}
    const unsigned short* getMismatchCyclesEnd() const {return mismatchCycles + mismatchCount;}
    void addMismatchCycle(const unsigned cycle)
    {
        ISAAC_ASSERT_MSG(cycle > 0, "Cycle numbers expected to be 1-based.");
        ISAAC_ASSERT_MSG(maxCycles_ >= cycle, "Cycle number is too high. Check maxCycles_.");
        mismatchCycles[mismatchCount] = cycle;
        ++mismatchCount;
    }


    std::string getCigarString() const
    {
        if (cigarBuffer && cigarLength)
        {
            return Cigar::toString(*cigarBuffer, cigarOffset, cigarLength);
        }
        else
        {
            return "";
        }
    }

    std::ostream &serializeCigar(std::ostream &os) const
    {
        if (cigarBuffer && cigarLength)
        {
            return Cigar::toStream(cigarBuffer->begin() + cigarOffset, cigarBuffer->begin() + cigarOffset + cigarLength, os);
        }
        else
        {
            return os;
        }
    }

    /**
     ** \brief Cluster associated to the fragment
     **
     ** Note, it is the responsibility of the calling code to ensure the
     ** validity of the cluster.
     **/
    const Cluster *cluster;
    /**
     ** \brief The cigarLength can be used to identify if a fragment has been
     ** successfully aligned
     **/
    bool isAligned() const {return 0 != cigarLength;}
    /**
     ** \brief Marks read as unaligned.
     **/
    void setUnaligned() {cigarBuffer = 0; cigarLength = 0; alignmentScore = -1U;}

    /**
     ** \brief Marks read as something that has no match position. This is different from setUnaligned
     **        as unaligned shadows still have a position of their orphan
     **/
    void setNoMatch()
        {setUnaligned(); contigId = reference::ReferencePosition::MAX_CONTIG_ID; position = 0;}
    bool isNoMatch() const {return reference::ReferencePosition::MAX_CONTIG_ID == contigId;}

    /**
     * \brief Notion of uniquely aligned in CASAVA means that a fragment was
     *        seen to have only one candidate alignment position. As this is
     *        highly dependent on the choice of seeds, alignment score based
     *        approximation should do.
     */
    bool isUniquelyAligned() const {return isAligned() && hasAlignmentScore() && 3 < getAlignmentScore();}

    unsigned getQuality() const
    {
        const std::vector<char> &quality = getRead().getForwardQuality();

        return std::accumulate(quality.begin(), quality.end(), 0U);
    }

    const std::vector<char> &getStrandSequence() const {return getRead().getStrandSequence(reverse);}
    const std::vector<char> &getStrandQuality() const {return getRead().getStrandQuality(reverse);}

    int getFirstSeedIndex() const {return firstSeedIndex;}

    bool hasAlignmentScore() const {return -1U != alignmentScore;}

    void incrementClipLeft(const unsigned short bases) {position += bases; if (reverse) {highClipped += bases;} else {lowClipped += bases;}}
    void incrementClipRight(const unsigned short bases) {if (reverse) {lowClipped += bases;} else {highClipped += bases;}}

    /// \return number of bases clipped on the left side of the fragment in respect to the reference
    unsigned short leftClipped() const {return reverse ? highClipped : lowClipped;}
    /// \return number of bases clipped on the right side of the fragment in respect to the reference
    unsigned short rightClipped() const {return reverse ? lowClipped : highClipped;}

    /// \return number of bases clipped on the left side of the fragment in respect to the reference
    unsigned short &leftClipped() {return reverse ? highClipped : lowClipped;}
    /// \return number of bases clipped on the right side of the fragment in respect to the reference
    unsigned short &rightClipped() {return reverse ? lowClipped : highClipped;}

    void resetAlignment(Cigar &buffer)
    {
        // this needs to be done before cigarOffset is reset;
        position = getUnclippedPosition();
        cigarOffset = buffer.size();
        cigarLength = 0;
        cigarBuffer = &buffer;
        observedLength = 0;
        std::fill(mismatchCycles, mismatchCycles + mismatchCount, 0);
        mismatchCount = 0;
        matchesInARow = 0;
        gapCount = 0;
        editDistance = 0;
        logProbability = 0.0;
        alignmentScore = -1U;
        smithWatermanScore = 0;
    }
    void resetClipping()
    {
        ISAAC_ASSERT_MSG(!isAligned(), "Alignment must be reset before clipping");
        lowClipped = 0;
        highClipped = 0;
    }

    void consolidate(const FragmentMetadata &that);

    bool isWellAnchored() const;

    unsigned getContigId() const {return contigId;}

    long getPosition() const {return position;}

    /// Id of the contig where the fragment is located
    unsigned contigId;
    /**
     ** \brief 0-based leftmost position of the fragment on the forward strand
     ** of the contig.
     **
     ** Even though the position can become negative while building the fragment
     ** (before calculating the cigar), the final position will be guaranteed to
     ** be positive (or 0). The final position is the position of the first
     ** aligned base ('M', '=' or 'X'). If the read extends outside the contig,
     ** this will be reflected by appropriate insertions or clipping operations
     ** in the cigar.
     **/
    long position;

    /// number of bases clipped from the lowest read cycle irrespective of alignment
    unsigned short lowClipped;
    /// number of bases clipped from the highest read cycle irrespective of alignment
    unsigned short highClipped;

    /**
     ** \brief observed length of the fragment on the contig
     **
     ** If there are no indels and no clipping, this would be the length of the
     ** read. With indels this is the read length minus the insertions plus the
     ** deletions (resp. to and from the reference). Clipped areas of the read
     ** are also subtracted from the length of the read.
     **/
    unsigned observedLength;
    /// 0-based index of the read in the list of ReadMetadata
    unsigned readIndex;
    /// Orientation of the read. False is forward, true is reverse
    bool reverse;
    /// Cigar offset in the associated cigar buffer (see FragmentBuilder)
    unsigned cigarOffset;
    /// Number of operations in the cigar
    unsigned cigarLength;
    /// Buffer containing the cigar data
    const std::vector<unsigned> *cigarBuffer;
    static const unsigned maxCycles_ = 1024;

    /// Number of mismatches in the alignment (can't be more than read length)
    unsigned mismatchCount;

    /// Longest stretch of matches
    unsigned matchesInARow;

    /// Number of short indels in the fragment
    unsigned gapCount;

    /// Edit distance from the alignment (including indels and ambiguous bases)
    unsigned editDistance;

    // array of cycle numbers containing mismatches (outside indels).
    unsigned short mismatchCycles[maxCycles_];
    /**
     ** \brief natural logarithm of the probability that the fragment is correct
     **
     ** This is intrinsic to the fragment and depends only of the quality of all
     ** matching base and and the quality of all mismatching bases (indel are
     ** counted as matching bases). See AlignmentQuality for the detail of the
     ** probabilities used in each case.
     **/
    double logProbability;

    /// The id of the seed that produced the alignment candidate. Valid only prior to consolidation.
    int firstSeedIndex;

    /// Count of seeds that mapped to highly repetitive locations in genome
    unsigned repeatSeedsCount;
    /// Count of seeds that mapped to that fragment that don't have neighbors in the reference within
    /// hamming distance of 4
    unsigned uniqueSeedCount;
    /// highest and lowest seed offsets for matches to the kmers having neighbors in the reference
    std::pair<unsigned, unsigned> nonUniqueSeedOffsets;
    /**
     ** \brief Alignment score in the global context of the reference
     **
     ** This depends on the all the pLogCorrect values for all the possible
     ** alignments for this read across the whole reference. It also takes into
     ** account the rest-of-genome correction. Value of -1U indicates unknown alignment score.
     **/
    unsigned alignmentScore;

    // Weighted sum of mismatch and gap penalties similar to what's used for Smith-Waterman alignment
    unsigned smithWatermanScore;

    /**
     ** \brief Comparison of FragmentMetadata by reference position
     **/
    bool operator<(const FragmentMetadata &f) const
    {
        ISAAC_ASSERT_MSG(cluster == f.cluster && readIndex == f.readIndex,
                         "Comparison makes sense only for metadata representing the same fragment " <<
                         *this << " vs " << f);
        return
            (this->contigId < f.contigId ||
                (this->contigId == f.contigId && (this->position < f.position ||
                    (this->position == f.position && (this->reverse < f.reverse ||
                        (reverse == f.reverse && observedLength < f.observedLength))))));
    }

    bool operator == (const FragmentMetadata &that) const
    {
        ISAAC_ASSERT_MSG(cluster == that.cluster && readIndex == that.readIndex,
                         "Comparison makes sense only for metadata representing the same fragment " <<
                         *this << " vs " << that);
        return position == that.position && contigId == that.contigId && reverse == that.reverse &&
            observedLength ==that.observedLength;
    }

    bool operator != (const FragmentMetadata &that) const
    {
        return !(*this == that);
    }
};

typedef std::vector<FragmentMetadata> FragmentMetadataList;

inline std::ostream &operator<<(std::ostream &os, const FragmentMetadata &f)
{
    os << "FragmentMetadata("
              << f.getCluster().getId() << "id "
              << f.contigId << ":"
              << f.position << ", "
              << f.observedLength << "bp, r"
              << f.readIndex << ", "
              << (f.reverse ? 'R' : 'F') << ", "
              << f.mismatchCount << "mm, "
              << f.matchesInARow << "mir, "
              << f.gapCount << "g, "
              << f.editDistance << "ed, ";
    return f.serializeCigar(os) << ", "
              << f.logProbability << "lp, "
              << f.repeatSeedsCount << "rs, "
              << f.uniqueSeedCount << "usc, "
              << f.alignmentScore << "sm, "
              << f.smithWatermanScore << "sws,"
              << f.isWellAnchored() << "wa" << ")";
}

inline void FragmentMetadata::consolidate(const FragmentMetadata &that)
{
    uniqueSeedCount += that.uniqueSeedCount;
    nonUniqueSeedOffsets.first = std::min(nonUniqueSeedOffsets.first, that.nonUniqueSeedOffsets.first);
    nonUniqueSeedOffsets.second = std::max(nonUniqueSeedOffsets.second, that.nonUniqueSeedOffsets.second);
}

inline bool FragmentMetadata::isWellAnchored() const
{
    return uniqueSeedCount ||
        (nonUniqueSeedOffsets.second > nonUniqueSeedOffsets.first &&
            // weak seeds must not overlap for well-anchored alignment
            (nonUniqueSeedOffsets.second - nonUniqueSeedOffsets.first) >= WEAK_SEED_LENGTH);
}



} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_METADATA_HH
