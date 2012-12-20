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
 ** \file SemialignedEndsClipper.cpp
 **
 ** \brief See SemialignedEndsClipper.hh
 ** 
 ** \author Roman Petrovski
 **/

#include "alignment/Alignment.hh"
#include "alignment/FragmentMetadata.hh"
#include "alignment/matchSelector/SemialignedEndsClipper.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

bool SemialignedEndsClipper::clipLeftSide(
    const std::vector<reference::Contig> &contigList,
    FragmentMetadata &fragmentMetadata)
{
    const Read &read = fragmentMetadata.getRead();

    std::vector<char>::const_iterator sequenceBegin = read.getStrandSequence(fragmentMetadata.reverse).begin();
    unsigned oldCigarOffset = fragmentMetadata.cigarOffset;
    unsigned oldCigarLength = fragmentMetadata.cigarLength;
    std::pair<unsigned, Cigar::OpCode> operation = Cigar::decode(fragmentMetadata.cigarBuffer->at(oldCigarOffset));
    unsigned softClippedBeginBases = 0;
    if (Cigar::SOFT_CLIP == operation.second)
    {
        if (2 > fragmentMetadata.cigarLength)
        {
            //when the adapter sequence
            // happens to be at the extremities of the read, the whole read gets clipped away.
            return false;
        }
        ++oldCigarOffset;
        --oldCigarLength;
        softClippedBeginBases = operation.first;
        sequenceBegin += operation.first;
        operation = Cigar::decode(fragmentMetadata.cigarBuffer->at(oldCigarOffset));
    }

    if (Cigar::ALIGN == operation.second)
    {
        unsigned mappedBeginBases = operation.first;
        std::vector<char>::const_iterator sequenceEnd = sequenceBegin + mappedBeginBases;

        const std::vector<char> &reference = contigList.at(fragmentMetadata.contigId).forward_;
        std::vector<char>::const_iterator referenceBegin = reference.begin() + fragmentMetadata.position;

        unsigned clipped = clipMismatches<CONSECUTIVE_MATCHES_MIN>(sequenceBegin, sequenceEnd,
                                                                   referenceBegin, reference.end(),
                                                                   &boost::cref<char>);

        if (clipped)
        {
            fragmentMetadata.cigarOffset = cigarBuffer_.size();
            fragmentMetadata.observedLength -= clipped;
            softClippedBeginBases += clipped;
            mappedBeginBases -= clipped;
            fragmentMetadata.position += clipped;

            cigarBuffer_.push_back(Cigar::encode(softClippedBeginBases, Cigar::SOFT_CLIP));
            cigarBuffer_.push_back(Cigar::encode(mappedBeginBases, Cigar::ALIGN));
            cigarBuffer_.insert(cigarBuffer_.end(),
                                fragmentMetadata.cigarBuffer->begin() + oldCigarOffset + 1,
                                fragmentMetadata.cigarBuffer->begin() + oldCigarOffset + oldCigarLength);
            fragmentMetadata.cigarBuffer = &cigarBuffer_;
            fragmentMetadata.cigarLength = cigarBuffer_.size() - fragmentMetadata.cigarOffset;

            return true;
        }
    }

    return false;
}

void SemialignedEndsClipper::clipRightSide(
    const std::vector<reference::Contig> &contigList,
    FragmentMetadata &fragmentMetadata)
{
    const Read &read = fragmentMetadata.getRead();

    std::reverse_iterator<std::vector<char>::const_iterator> sequenceRBegin(
        read.getStrandSequence(fragmentMetadata.reverse).end());
    unsigned oldCigarOffset = fragmentMetadata.cigarOffset;
    unsigned oldCigarLength = fragmentMetadata.cigarLength;
    std::pair<unsigned, Cigar::OpCode> operation = Cigar::decode(
        fragmentMetadata.cigarBuffer->at(oldCigarOffset + oldCigarLength - 1));
    unsigned softClippedEndBases = 0;
    if (Cigar::SOFT_CLIP == operation.second)
    {
        if (2 > fragmentMetadata.cigarLength)
        {
            //when the adapter sequence
            // happens to be at the extremities of the read, the whole read gets clipped away.
            return;
        }

        --oldCigarLength;
        softClippedEndBases = operation.first;
        sequenceRBegin += operation.first;
        operation = Cigar::decode(fragmentMetadata.cigarBuffer->at(oldCigarOffset + oldCigarLength - 1));
    }

    if (Cigar::ALIGN == operation.second)
    {
        unsigned mappedEndBases = operation.first;
        std::reverse_iterator<std::vector<char>::const_iterator> sequenceREnd = sequenceRBegin + mappedEndBases;

        const std::vector<char> &reference = contigList.at(fragmentMetadata.contigId).forward_;
        std::reverse_iterator<std::vector<char>::const_iterator> referenceRBegin(reference.begin() + fragmentMetadata.position +
            fragmentMetadata.getObservedLength());
        std::reverse_iterator<std::vector<char>::const_iterator> referenceREnd(reference.begin());

        unsigned clipped = clipMismatches<CONSECUTIVE_MATCHES_MIN>(sequenceRBegin, sequenceREnd,
                                                                   referenceRBegin, referenceREnd,
                                                                   &boost::cref<char>);

        if (clipped)
        {
            fragmentMetadata.cigarOffset = cigarBuffer_.size();
            fragmentMetadata.observedLength -= clipped;
            softClippedEndBases += clipped;
            mappedEndBases -= clipped;
            cigarBuffer_.insert(cigarBuffer_.end(),
                                fragmentMetadata.cigarBuffer->begin() + oldCigarOffset,
                                fragmentMetadata.cigarBuffer->begin() + oldCigarOffset + oldCigarLength - 1);
            cigarBuffer_.push_back(Cigar::encode(mappedEndBases, Cigar::ALIGN));
            cigarBuffer_.push_back(Cigar::encode(softClippedEndBases, Cigar::SOFT_CLIP));
            fragmentMetadata.cigarBuffer = &cigarBuffer_;
            fragmentMetadata.cigarLength = cigarBuffer_.size() - fragmentMetadata.cigarOffset;
        }
    }
}

/**
 * \return true, if the clipping changed the alignment position
 */
bool SemialignedEndsClipper::clip(
    const std::vector<reference::Contig> &contigList,
    FragmentMetadata &fragmentMetadata)
{
    if (!fragmentMetadata.isAligned())
    {
        return false;
    }

    bool ret = clipLeftSide(contigList, fragmentMetadata);
    clipRightSide(contigList, fragmentMetadata);

    return ret;
}

void SemialignedEndsClipper::clip(
    const std::vector<reference::Contig> &contigList,
    BamTemplate &bamTemplate)
{
    for (unsigned k = 0; k < bamTemplate.getFragmentCount(); ++k)
    {
        FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(k);
        if (clip(contigList, fragment))
        {
            if (2 == bamTemplate.getFragmentCount())
            {
                FragmentMetadata &mate = bamTemplate.getMateFragmentMetadata(fragment);
                if (!mate.isAligned())
                {
                    // For shadow mates, the position needs to be updated in case the soft clipping
                    // changed the position of the singleton.
                    mate.position = fragment.position;
                    // paired template mate is unaligned, no point to continue;
                    break;
                }
            }
        }
    }
}

} // namespace matchSelector
} // namespace alignment
} // namespace isaac
