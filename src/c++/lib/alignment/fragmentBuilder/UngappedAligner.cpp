/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** BSD 2-Clause License
 **
 ** You should have received a copy of the BSD 2-Clause License
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** \file UngappedAligner.cpp
 **
 ** \brief See UngappedAligner.hh
 ** 
 ** \author Roman Petrovski
 **/
#include "alignment/fragmentBuilder/UngappedAligner.hh"
#include "alignment/Alignment.hh"

namespace isaac
{
namespace alignment
{
namespace fragmentBuilder
{

UngappedAligner::UngappedAligner(
    const int gapMatchScore,
    const int gapMismatchScore,
    const int gapOpenScore,
    const int gapExtendScore,
    const int minGapExtendScore)
    : AlignerBase(gapMatchScore, gapMismatchScore, gapOpenScore, gapExtendScore, minGapExtendScore)
{
}

unsigned UngappedAligner::alignUngapped(
    FragmentMetadata &fragmentMetadata,
    Cigar &cigarBuffer,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::FragmentSequencingAdapterClipper &adapterClipper,
    const reference::Contig &contig) const
{
    const unsigned cigarOffset = cigarBuffer.size();

    fragmentMetadata.resetAlignment(cigarBuffer);
    fragmentMetadata.resetClipping();

    const Read &read = fragmentMetadata.getRead();
    const bool reverse = fragmentMetadata.reverse;
    const std::vector<char> &sequence = read.getStrandSequence(reverse);
    const std::vector<char> &reference = contig.forward_;

    std::vector<char>::const_iterator sequenceBegin = sequence.begin();
    std::vector<char>::const_iterator sequenceEnd = sequence.end();

    adapterClipper.clip(contig, fragmentMetadata, sequenceBegin, sequenceEnd);
    clipReadMasking(read, fragmentMetadata, sequenceBegin, sequenceEnd);

    clipReference(reference.size(), fragmentMetadata, sequenceBegin, sequenceEnd);

    const unsigned firstMappedBaseOffset = std::distance(sequence.begin(), sequenceBegin);
    if (firstMappedBaseOffset)
    {
        cigarBuffer.addOperation(firstMappedBaseOffset, Cigar::SOFT_CLIP);
    }

    const unsigned mappedBases = std::distance(sequenceBegin, sequenceEnd);
    if (mappedBases)
    {
        const Cigar::OpCode opCode = Cigar::ALIGN;
        cigarBuffer.addOperation(mappedBases, opCode);
    }

    const unsigned clipEndBases = std::distance(sequenceEnd, sequence.end());
    if (clipEndBases)
    {
        cigarBuffer.addOperation(clipEndBases, Cigar::SOFT_CLIP);
    }

    const unsigned ret = updateFragmentCigar(readMetadataList, reference, fragmentMetadata,
                                             fragmentMetadata.position, cigarBuffer, cigarOffset);

    if (!ret)
    {
        fragmentMetadata.setUnaligned();
    }

    return ret;
}

} // namespace fragmentBuilder
} // namespace alignment
} // namespace isaac
