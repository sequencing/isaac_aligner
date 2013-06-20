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
 ** \file SimpleIndelAligner.hh
 **
 ** \brief Uses seed match discrepancies to detect indels in the reads that span single event
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_SIMPLE_INDEL_ALIGNER_HH
#define iSAAC_ALIGNMENT_FRAGMENT_BUILDER_SIMPLE_INDEL_ALIGNER_HH


#include "alignment/fragmentBuilder/AlignerBase.hh"

namespace isaac
{
namespace alignment
{

namespace fragmentBuilder
{

class SimpleIndelAligner: public AlignerBase
{
    static const unsigned GAP_FLANK_BASES = 32;
    static const unsigned GAP_FLANK_MISMATCHES_MAX = 8;

    const unsigned semialignedGapLimit_;
public:
    SimpleIndelAligner(
        // defaults for unit tests
        const int gapMatchScore,
        const int gapMismatchScore,
        const int gapOpenScore,
        const int gapExtendScore,
        const int minGapExtendScore,
        const unsigned semialignedGapLimit);

    void alignSimpleIndels(
        Cigar &cigarBuffer,
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        const SeedMetadataList &seedMetadataList,
        FragmentMetadataList &fragmentList) const;

private:
    void alignSimpleDeletion(
        Cigar &cigarBuffer,
        FragmentMetadata &headFragment,
        const unsigned headSeedOffset,
        FragmentMetadata &tailFragment,
        const unsigned tailSeedOffset,
        const unsigned tailSeedLength,
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList) const;

    void alignSimpleInsertion(
        Cigar &cigarBuffer,
        FragmentMetadata &headAlignment,
        const unsigned headSeedOffset,
        const unsigned headSeedLength,
        FragmentMetadata &tailAlignment,
        const unsigned tailSeedOffset,
        const unsigned tailSeedLength,
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList) const;

};

} // namespace fragmentBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_SIMPLE_INDEL_ALIGNER_HH
