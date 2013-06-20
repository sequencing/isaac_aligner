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
 ** \file UngappedAligner.hh
 **
 ** \brief Aligns read to the reference without gaps
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_UNGAPPED_ALIGNER_HH
#define iSAAC_ALIGNMENT_FRAGMENT_BUILDER_UNGAPPED_ALIGNER_HH

#include "alignment/fragmentBuilder/AlignerBase.hh"

namespace isaac
{
namespace alignment
{

namespace fragmentBuilder
{

class Cluster;

/**
 ** \brief Utility component creating and scoring all Fragment instances from a
 ** list Seed Matches for a single Cluster (each Read independently).
 **/
class UngappedAligner: public AlignerBase
{
public:
    UngappedAligner(
        // defaults for unit tests
        const int gapMatchScore,
        const int gapMismatchScore,
        const int gapOpenScore,
        const int gapExtendScore,
        const int minGapExtendScore);

    /**
     ** \brief Calculate the ungapped alignment of a fragment
     **
     ** \return number of mapped matches. If 0, read is marked as unaligned
     **/
    unsigned alignUngapped(
        FragmentMetadata &fragmentMetadata,
        Cigar &cigarBuffer,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::FragmentSequencingAdapterClipper  &adapterClipper,
        const reference::Contig &contig) const;
};

} // namespace fragmentBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_UNGAPPED_ALIGNER_HH
