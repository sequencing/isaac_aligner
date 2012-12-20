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
 ** \file ShadowAligner.hh
 **
 ** \brief Aligns shadows
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_SHADOW_ALIGNER_HH
#define iSAAC_ALIGNMENT_SHADOW_ALIGNER_HH

#include <vector>
#include <boost/noncopyable.hpp>

#include "reference/Contig.hh"
#include "alignment/BandedSmithWaterman.hh"
#include "alignment/Cigar.hh"
#include "alignment/Cluster.hh"
#include "alignment/FragmentBuilder.hh"
#include "alignment/TemplateLengthStatistics.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Utility component aligning shadows.
 **
 ** The intended use case is for the TemplateBuilder to delegate the alignment
 ** of shadow reads (or poorly aligned mates) to this specialized component.
 **
 **/
class ShadowAligner: boost::noncopyable
{
public:
    /**
     ** \brief Construct a ShadowAligner for a reference genome and a given
     ** set of reads.
     **/
    ShadowAligner(const flowcell::FlowcellLayoutList &flowcellLayoutList,
                  const unsigned gappedMismatchesMax,
                  const FragmentBuilder &fragmentBuilder);
    /**
     ** \brief Helper method to align the shadow of an orphan.
     **
     ** Uses the position, read index, and orientation of the orphan to infers
     ** the orientation and range of positions of the shadow (according to the
     ** templateLengthStatistics_).
     **
     ** Hashes the all k-mers of length orphanKmerLength_ from the orphan
     ** sequence into the orphaKmerPositions_.
     **
     ** Find all the possible alignment positions of the orphan by comparing the
     ** k-mers of length orphanKmerLength_ from the region of the reference
     ** where the orphan could align to the hashed k-mers from the orphan
     ** sequence.
     **
     ** Align the orphan to the candidate positions.
     **
     ** \param orphan[in] the correctly aligned orphan
     **
     ** \param shadow[out] the unaligned shadow
     **/
    bool rescueShadow(
        const std::vector<reference::Contig> &contigList,
        const FragmentMetadata &orphan,
        std::vector<FragmentMetadata> &shadowList,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const TemplateLengthStatistics &templateLengthStatistics);
    const Cigar &getCigarBuffer() const {return shadowCigarBuffer_;}
private:
    static const unsigned unreasonablyHighDifferenceBetweenMaxAndMinInsertSizePlusFlanks_ = 10000;

    /// Length of the k-mers used to rescue shadows and misaligned reads
    static const unsigned shadowKmerLength_ = 7;
    static const unsigned shadowKmerCount_ = (1 << (2 * shadowKmerLength_));
    const unsigned gappedMismatchesMax_;
    const FragmentBuilder &fragmentBuilder_;
    /**
     ** \brief Cached storage for the position of the k-mers in the shadow
     **
     ** Note that this is a really fast and cheap but imperfect to rescue
     ** shadows or mis-aligned reads. The index used to access elements in the
     ** vector is made from a k-mer of length shadowKmerLength_ (the vector has
     ** 4 ^ shadowKmerLength_ positions). The values in the table at position i
     ** is the first position in the read where the k-mer was found (-1 if not
     ** found). Repeats are recorded only once in the table. This allows to
     ** identify extremely quickly if a k-mer in the reference belongs to the
     ** read. The shadowKmerLength_ should stay small enough to ensure that the
     ** table stays in the L1 cache.
     **/
    std::vector<short> shadowKmerPositions_;
    Cigar shadowCigarBuffer_;
    /// Hash all the k-mers of length shadowKmerLength_ into shadowKmerPositions_
    unsigned hashShadowKmers(const std::vector<char> &sequence);
    /**
     ** \brief Cached storage for the candidate start positions of the shadow
     **
     ** Note: all these positions are relative to the beginning of the region
     ** where the shadow is expected to be (as opposed to the position relative
     ** to the beginning of the reference).
     **/
    std::vector<long> shadowCandidatePositions_;
    /// Find all candidate positions for a shadow sequence on a given reference interval
    void findShadowCandidatePositions(
        const std::vector<char>::const_iterator referenceBegin,
        const std::vector<char>::const_iterator referenceEnd,
        const std::vector<char> &shadowSequence);
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SHADOW_ALIGNER_HH
