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
 ** \file SemialignedEndsClipper.cpp
 **
 ** \brief See SemialignedEndsClipper.hh
 ** 
 ** \author Roman Petrovski
 **/

#include "alignment/Alignment.hh"
#include "common/Debug.hh"

#include "SemialignedEndsClipper.hh"

namespace isaac
{
namespace build
{

using alignment::Cigar;

/**
 * \brief clips mismatches on the left if this does not move index.pos_ to binEndPos or beyond
 */
bool SemialignedEndsClipper::clipLeftSide(
    const std::vector<reference::Contig> &contigList,
    const reference::ReferencePosition binEndPos,
    PackedFragmentBuffer::Index &index,
    io::FragmentAccessor &fragment)
{
    const unsigned char *sequenceBegin = fragment.basesBegin();

    const uint32_t *oldCigarBegin = index.cigarBegin_;
    std::pair<unsigned, alignment::Cigar::OpCode> operation = alignment::Cigar::decode(*oldCigarBegin);
    unsigned softClippedBeginBases = 0;
    if (alignment::Cigar::SOFT_CLIP == operation.second)
    {
        ++oldCigarBegin;
        softClippedBeginBases = operation.first;
        sequenceBegin += operation.first;
        operation = alignment::Cigar::decode(*oldCigarBegin);
    }

    if (alignment::Cigar::ALIGN == operation.second)
    {
        unsigned mappedBeginBases = operation.first;
        const unsigned char * sequenceEnd = sequenceBegin + mappedBeginBases;

        const std::vector<char> &reference = contigList.at(index.pos_.getContigId()).forward_;
        std::vector<char>::const_iterator referenceBegin = reference.begin() + index.pos_.getPosition();

        std::pair<unsigned, unsigned> clipped = alignment::clipMismatches<CONSECUTIVE_MATCHES_MIN>(sequenceBegin, sequenceEnd,
                                                                              referenceBegin, reference.end(),
                                                                              &oligo::getUppercaseBaseFromBcl);

        if (clipped.first && index.pos_ + clipped.first < binEndPos)
        {
            softClippedBeginBases += clipped.first;
            mappedBeginBases -= clipped.first;
            index.pos_ += clipped.first;
            fragment.fStrandPosition_ += clipped.first;
            fragment.observedLength_ -= clipped.first;
            fragment.editDistance_ -= clipped.second;

            size_t before = cigarBuffer_.size();
            cigarBuffer_.push_back(Cigar::encode(softClippedBeginBases, Cigar::SOFT_CLIP));
            cigarBuffer_.push_back(Cigar::encode(mappedBeginBases, Cigar::ALIGN));
            cigarBuffer_.insert(cigarBuffer_.end(),
                                oldCigarBegin + 1,
                                index.cigarEnd_);
            index.cigarBegin_ = &cigarBuffer_.at(before);
            index.cigarEnd_ = &cigarBuffer_.back() + 1;
            return true;
        }
    }
    return false;
}

bool SemialignedEndsClipper::clipRightSide(
    const std::vector<reference::Contig> &contigList,
    PackedFragmentBuffer::Index &index,
    io::FragmentAccessor &fragment)
{
    const unsigned char *sequenceEnd = fragment.basesEnd();

    std::reverse_iterator<const unsigned char *> sequenceRBegin(sequenceEnd);

    const uint32_t *oldCigarEnd = index.cigarEnd_;
    std::pair<unsigned, alignment::Cigar::OpCode> operation = alignment::Cigar::decode(*(oldCigarEnd - 1));
    unsigned softClippedEndBases = 0;
    if (Cigar::SOFT_CLIP == operation.second)
    {
        --oldCigarEnd;
        softClippedEndBases = operation.first;
        sequenceRBegin += operation.first;
        operation = alignment::Cigar::decode(*(oldCigarEnd - 1));
    }

    if (alignment::Cigar::ALIGN == operation.second)
    {
        unsigned mappedEndBases = operation.first;
        std::reverse_iterator<const unsigned char *> sequenceREnd = sequenceRBegin + mappedEndBases;

        const std::vector<char> &reference = contigList.at(index.pos_.getContigId()).forward_;
        std::reverse_iterator<std::vector<char>::const_iterator> referenceRBegin(reference.begin() + index.pos_.getPosition() +
            fragment.observedLength_);
        std::reverse_iterator<std::vector<char>::const_iterator> referenceREnd(reference.begin());

        std::pair<unsigned, unsigned> clipped = alignment::clipMismatches<CONSECUTIVE_MATCHES_MIN>(sequenceRBegin, sequenceREnd,
                                                                              referenceRBegin, referenceREnd,
                                                                              &oligo::getUppercaseBaseFromBcl);

        if (clipped.first)
        {
            softClippedEndBases += clipped.first;
            mappedEndBases -= clipped.first;
            fragment.observedLength_ -= clipped.first;
            fragment.editDistance_ -= clipped.second;

            size_t before = cigarBuffer_.size();
            cigarBuffer_.insert(cigarBuffer_.end(), index.cigarBegin_, oldCigarEnd - 1);
            cigarBuffer_.push_back(Cigar::encode(mappedEndBases, Cigar::ALIGN));
            cigarBuffer_.push_back(Cigar::encode(softClippedEndBases, Cigar::SOFT_CLIP));
            index.cigarBegin_ = &cigarBuffer_.at(before);
            index.cigarEnd_ = &cigarBuffer_.back() + 1;
            return true;
        }
    }
    return false;
}

void SemialignedEndsClipper::clip(
    const std::vector<reference::Contig> &contigs,
    const reference::ReferencePosition binEndPos,
    PackedFragmentBuffer::Index &index,
    io::FragmentAccessor &fragment)
{
    ISAAC_ASSERT_MSG(fragment.isAligned(), "Unexpected unaligned fragment from gap realigner");

    const bool leftClipped = clipLeftSide(contigs, binEndPos, index, fragment);
    const bool rightClipped = clipRightSide(contigs, index, fragment);
    if (leftClipped || rightClipped)
    {
        ISAAC_THREAD_CERR_DEV_TRACE(" SemialignedEndsClipper::clip: " << fragment);
    }
}

} // namespace build
} // namespace isaac
