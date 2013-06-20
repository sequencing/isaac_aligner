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
 ** \file GappedAligner.hh
 **
 ** \brief Uses banded Smith-Waterman algorithm to align fragment
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_GAPPED_ALIGNER_HH
#define iSAAC_ALIGNMENT_FRAGMENT_BUILDER_GAPPED_ALIGNER_HH

#include "alignment/fragmentBuilder/AlignerBase.hh"
#include "alignment/BandedSmithWaterman.hh"

namespace isaac
{
namespace alignment
{

namespace fragmentBuilder
{

class GappedAligner: public AlignerBase
{
public:
    GappedAligner(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const bool avoidSmithWaterman,
        const int gapMatchScore,
        const int gapMismatchScore,
        const int gapOpenScore,
        const int gapExtendScore,
        const int minGapExtendScore);

    /**
     ** \brief Calculate the gapped alignment of a fragment
     **/
    unsigned alignGapped(
        FragmentMetadata &fragmentMetadata,
        Cigar &cigarBuffer,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::FragmentSequencingAdapterClipper &adapterClipper,
        const reference::Contig &contig);

protected:
    const bool avoidSmithWaterman_;
    BandedSmithWaterman bandedSmithWaterman_;
    static const unsigned HASH_KMER_LENGTH = 7;
    static const unsigned QUERY_LENGTH_MAX = 65536;

    common::FiniteCapacityVector<unsigned, 2> hashedQueryTile_;
    common::FiniteCapacityVector<unsigned, 2> hashedQueryCluster_;
    common::FiniteCapacityVector<unsigned, 2> hashedQueryReadIndex_;

//    std::vector<char>::const_iterator hashedQueryBegin_;
//    std::vector<char>::const_iterator hashedQueryEnd_;
    // initialize all k-mers to the magic value -1 (NOT_FOUND)
    common::FiniteCapacityVector<unsigned short, oligo::MaxKmer<HASH_KMER_LENGTH, unsigned short>::value + 1> queryKmerOffsets_;
    static const unsigned short UNINITIALIZED_OFFSET_MAGIC = static_cast<unsigned short>(-1);
    static const unsigned short REPEAT_OFFSET_MAGIC = static_cast<unsigned short>(-2);
    // count of hits required to assume that part of the sequence will anchor at a position.
    // the smaller the number, the more reads will go into smith-waterman. The higher the number
    // the more likely is that some good gapped alignments will not be attempted
    static const unsigned char SUFFICIENT_NUMBER_OF_HITS = 8;



    bool makesSenseToGapAlign(
        const unsigned tile, const unsigned cluster, const unsigned read, const bool reverse,
        const std::vector<char>::const_iterator queryBegin,
        const std::vector<char>::const_iterator queryEnd,
        const std::vector<char>::const_iterator databaseBegin,
        const std::vector<char>::const_iterator databaseEnd);
};

} // namespace fragmentBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_GAPPED_ALIGNER_HH
