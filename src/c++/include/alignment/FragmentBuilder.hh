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
 ** \file FragmentBuilder.hh
 **
 ** \brief Utility classes for Fragment building and management for several reads
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_HH
#define iSAAC_ALIGNMENT_FRAGMENT_BUILDER_HH

#include "alignment/fragmentBuilder/GappedAligner.hh"
#include "alignment/fragmentBuilder/UngappedAligner.hh"
#include "alignment/fragmentBuilder/SimpleIndelAligner.hh"
#include "alignment/BandedSmithWaterman.hh"
#include "alignment/Cigar.hh"
#include "alignment/FragmentMetadata.hh"
#include "alignment/Match.hh"
#include "alignment/SeedMetadata.hh"
#include "alignment/matchSelector/FragmentSequencingAdapterClipper.hh"
#include "reference/Contig.hh"
#include "flowcell/ReadMetadata.hh"

namespace isaac
{
namespace alignment
{

class Cluster;

/**
 ** \brief Utility component creating and scoring all Fragment instances from a
 ** list Seed Matches for a single Cluster (each Read independently).
 **/
class FragmentBuilder: public boost::noncopyable
{
public:
    FragmentBuilder(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned repeatThreshold,
        const unsigned maxSeedsPerRead,
        const unsigned gappedMismatchesMax,
        const bool avoidSmithWaterman,
        const int gapMatchScore,
        const int gapMismatchScore,
        const int gapOpenScore,
        const int gapExtendScore,
        const int minGapExtendScore,
        const unsigned semialignedGapLimit);

    bool build(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        const SeedMetadataList &seedMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<Match>::const_iterator matchBegin,
        const std::vector<Match>::const_iterator matchEnd,
        const Cluster &cluster,
        bool withGaps);
    const std::vector<std::vector<FragmentMetadata> > &getFragments() const {return fragments_;}
    const std::vector<unsigned> &getCigarBuffer() const {return cigarBuffer_;}

    struct SequencingAdapterRange
    {
        SequencingAdapterRange() : defined_(false), empty_(true){}
        bool defined_;
        bool empty_;
        std::vector<char>::const_iterator adapterRangeBegin_;
        std::vector<char>::const_iterator adapterRangeEnd_;
    };

private:
    static const unsigned readsMax_ = 2;
    const unsigned repeatThreshold_;
    const unsigned semialignedGapLimit_;
    const unsigned gappedMismatchesMax_;

    /**
     * \brief flag per seed indicating whether the seed matches are ignored due to
     *        a high repeat match
     */
    std::vector<unsigned> seedMatchCounts_;
    unsigned repeatSeedsCount_;
    /**
     ** \brief All FragmentMetadata for all reads
     **
     ** fragments_[i] is the list of fragments for read i.
     **/
    std::vector<std::vector<FragmentMetadata> > fragments_;
    /**
     ** \brief Buffer for all the CIGAR
     **
     ** The buffer store all the CIGAR for all the FragmentMetadata. This is
     ** done to avoid allocating memory every time a new CIGAR is created.
     **/
    Cigar cigarBuffer_;

    fragmentBuilder::UngappedAligner ungappedAligner_;
    fragmentBuilder::GappedAligner gappedAligner_;
    fragmentBuilder::SimpleIndelAligner simpleIndelAligner_;

    /// clear all the buffers
    void clear();
    /**
     ** \brief add a match, either by creating a new instance of
     ** FragmentMetadata or by updating an existing one
     **
     ** Initializes the FragmentMetadata in the list for the corresponding
     ** readIndex with contigId, orientation (reverse flag) and
     ** position. The fragment is initially located at the leftmost position of
     ** the read on the forward strand of the contig. This means that the
     ** position can be negative.
     **
     ** Note: spurious FragmentMetadata creation is avoided by checking if the
     ** last FragmentMetadata created for the read has same contigId, position
     ** and orientation.
     **/
    void addMatch(
        const flowcell::ReadMetadataList &readMetadataList,
        const SeedMetadataList &seedMetadataList, const Match &match,
        const Cluster &cluster);
    /// Calculate the alignment for all fragments identified so far
    void alignFragments(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        const SeedMetadataList &seedMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        bool withGaps);

    /**
     ** \brief Position of the leftmost base of a read on the forward strand,
     ** given a seed, its position and orientation.
     **
     ** For forward matches, the offset to apply is simply the seed offset
     ** indicated in the SeedMetadata. For reverse matches, the offset to apply
     ** is the remaining length of the read after the end of the seed.
     **/
    long getReadPosition(
        const flowcell::ReadMetadataList &readMetadataList,
        const SeedMetadata &seedMetadata,
        const long seedPosition,
        const bool reverse) const;

    /// consolidate fragments with same reference position and orientation for a single read
    static void consolidateDuplicateFragments(FragmentMetadataList &fragmentList, const bool removeUnaligned);
    void removeRepeatSeedAlignments(FragmentMetadataList &fragmentList);
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_HH
