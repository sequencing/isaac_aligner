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
 ** \file OverlappingEndsClipper.cpp
 **
 ** \brief See OverlappingEndsClipper.hh
 ** 
 ** \author Roman Petrovski
 **/

#include "alignment/Alignment.hh"
#include "alignment/FragmentMetadata.hh"
#include "alignment/matchSelector/OverlappingEndsClipper.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{


/**
 * \return true, if the clipping changed the alignment position
 */
int OverlappingEndsClipper::diffBaseQualities(
    std::vector<char>::const_iterator left,
    std::vector<char>::const_iterator right,
    unsigned length)
{
    const int ret = std::inner_product( left, left + length, right, 0, std::plus<int>(), std::minus<int>());
    return ret;
}

void OverlappingEndsClipper::clip(
    const std::vector<reference::Contig> &contigList,
    BamTemplate &bamTemplate)
{
    if (2 != bamTemplate.getFragmentCount())
    {
        return;
    }

    FragmentMetadata &r1 = bamTemplate.getFragmentMetadata(0);
    FragmentMetadata &r2 = bamTemplate.getFragmentMetadata(1);

    if (!r1.isAligned() || !r2.isAligned())
    {
        return;
    }

    if (r1.contigId != r1.contigId)
    {
        // ignore chimeric pairs
        return;
    }

    if (r1.isReverse() == r2.isReverse())
    {
        return;
    }

    FragmentMetadata &left = r1.position < r2.position ? r1 : r2;
    FragmentMetadata &right = r1.position <= r2.position ? r2 : r1;

    if (left.isReverse())
    {
        // if leftmost one is not facing backwards, it's either a mate-pair or we've sequenced the adapter.
        // Deal with it somewhere else...
        return;
    }

    const long overlapLength = left.position + left.getObservedLength() - right.position;
    if (0 >= overlapLength)
    {
        // no overlap
        return;
    }

    // get the overlapping bit of the left read
    unsigned leftEndSoftClip = 0;
    unsigned leftEndOffset = left.getReadLength();
    Cigar::const_iterator leftLastCigarAlignOpIt = left.cigarEnd() - 1;
    Cigar::Component leftLastCigarAlignOp = Cigar::decode(*leftLastCigarAlignOpIt);
    if (Cigar::SOFT_CLIP == leftLastCigarAlignOp.second)
    {
        ISAAC_ASSERT_MSG(left.cigarBegin() != leftLastCigarAlignOpIt, "Fully soft-clipped reads are not allowed");
        leftEndOffset -= leftLastCigarAlignOp.first;
        leftEndSoftClip = leftLastCigarAlignOp.first;
        --leftLastCigarAlignOpIt;
        leftLastCigarAlignOp = Cigar::decode(*leftLastCigarAlignOpIt);
    }

    ISAAC_ASSERT_MSG(Cigar::ALIGN == leftLastCigarAlignOp.second,
                     "Apart from soft-clipping, CIGAR must end with align operations." << left);

    if (overlapLength >= leftLastCigarAlignOp.first)
    {
        // overlap contains or borders an indel or read will get fully soft-clipped, don't mess with those.
        return;
    }

    // get the overlapping bit of the right read
    unsigned rightStartOffset = 0;
    Cigar::const_iterator rightFirstCigarAlignOpIt = right.cigarBegin();
    Cigar::Component rightFirstCigarAlignOp = Cigar::decode(*rightFirstCigarAlignOpIt);
    if (Cigar::SOFT_CLIP == rightFirstCigarAlignOp.second)
    {
        rightStartOffset += rightFirstCigarAlignOp.first;
        ++rightFirstCigarAlignOpIt;
        ISAAC_ASSERT_MSG(right.cigarEnd() != rightFirstCigarAlignOpIt, "Fully soft-clipped reads are not allowed");
        rightFirstCigarAlignOp = Cigar::decode(*rightFirstCigarAlignOpIt);
    }

    ISAAC_ASSERT_MSG(Cigar::ALIGN == rightFirstCigarAlignOp.second,
                     "Apart from soft-clipping, CIGAR must begin with align operations." << right);

    if (overlapLength >= rightFirstCigarAlignOp.first)
    {
        // overlap contains or borders an indel or read will get fully soft-clipped, don't mess with those.
        return;
    }

    // find which of the overlapping ends is better
    if (0 < diffBaseQualities(
        left.getRead().getForwardQuality().begin() + leftEndOffset - overlapLength,
        right.getRead().getReverseQuality().begin() + rightStartOffset, overlapLength))
    {
        // left one is better, clip right
        const std::vector<char>::const_iterator reference =
            contigList.at(right.contigId).forward_.begin() + right.position;

        // right one is better
        Cigar::const_iterator rightCigarEndIt = right.cigarEnd();
        right.cigarOffset = cigarBuffer_.size();
        cigarBuffer_.push_back(Cigar::encode(rightStartOffset + overlapLength, Cigar::SOFT_CLIP));
        cigarBuffer_.push_back(Cigar::encode(rightFirstCigarAlignOp.first - overlapLength, Cigar::ALIGN));
        cigarBuffer_.insert(cigarBuffer_.end(), rightFirstCigarAlignOpIt + 1, rightCigarEndIt);
        right.incrementClipLeft(overlapLength);

        right.observedLength -= overlapLength;
        right.editDistance -= std::inner_product(
            right.getRead().getReverseSequence().begin() + rightStartOffset,
            right.getRead().getReverseSequence().begin() + rightStartOffset + overlapLength, reference,
            0, std::plus<int>(), std::not_equal_to<char>());
        right.cigarBuffer = &cigarBuffer_;
        right.cigarLength = cigarBuffer_.size() - right.cigarOffset;
    }
    else
    {
        // right one is better, clip left
        const std::vector<char>::const_iterator reference =
            contigList.at(left.contigId).forward_.begin() +
            left.position + left.getObservedLength() - leftEndSoftClip - overlapLength;

        const Cigar::const_iterator leftCigarBeginIt = left.cigarBegin();
        left.cigarOffset = cigarBuffer_.size();
        cigarBuffer_.insert(cigarBuffer_.end(), leftCigarBeginIt, leftLastCigarAlignOpIt);
        cigarBuffer_.push_back(Cigar::encode(leftLastCigarAlignOp.first - overlapLength, Cigar::ALIGN));
        cigarBuffer_.push_back(Cigar::encode(leftEndSoftClip + overlapLength, Cigar::SOFT_CLIP));
        left.incrementClipRight(overlapLength);

        left.observedLength -= overlapLength;
        left.editDistance -= std::inner_product(
            left.getRead().getForwardSequence().begin() + leftEndOffset - overlapLength,
            left.getRead().getForwardSequence().begin() + leftEndOffset, reference,
            0, std::plus<int>(), std::not_equal_to<char>());
        left.cigarBuffer = &cigarBuffer_;
        left.cigarLength = cigarBuffer_.size() - left.cigarOffset;

    }
}

} // namespace matchSelector
} // namespace alignment
} // namespace isaac
