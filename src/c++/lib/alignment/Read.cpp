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
 ** \file Read.cpp
 **
 ** Component containing the data associated to a cluster: sequence and quality
 ** strings for all the reads in the cluster.
 **
 ** \author Come Raczy
 **/

#include <cassert>
#include <boost/foreach.hpp>

#include "oligo/Nucleotides.hh"
#include "alignment/Read.hh"

namespace isaac
{
namespace alignment
{

void Read::decodeBcl(
    std::vector<char>::const_iterator bclBegin,
    std::vector<char>::const_iterator bclEnd,
    const unsigned index)
{
    ISAAC_ASSERT_MSG(index_ == index, "Unexpected initialization of one read with another");
    assert(bclEnd > bclBegin);
    const size_t readLength = bclEnd - bclBegin;
    ISAAC_ASSERT_MSG(forwardSequence_.capacity() >= readLength, "Buffers expected to be preallocated");
    ISAAC_ASSERT_MSG(reverseSequence_.capacity() >= readLength, "Buffers expected to be preallocated");
    ISAAC_ASSERT_MSG(forwardQuality_.capacity() >= readLength, "Buffers expected to be preallocated");
    ISAAC_ASSERT_MSG(reverseQuality_.capacity() >= readLength, "Buffers expected to be preallocated");
//    sequence_.reserve(readLength);
    forwardSequence_.clear();
    reverseSequence_.clear();
//    quality_.reserve(readLength);
    forwardQuality_.clear();
    reverseQuality_.clear();
    /*beginCyclesMasked_ = 0L;*/
    endCyclesMasked_ = 0L;

    for (std::vector<char>::const_iterator bcl = bclBegin; bclEnd > bcl; ++bcl)
    {

        if (!oligo::isBclN(*bcl))
        {
            forwardSequence_.push_back(oligo::getBase((*bcl) & 3, true));
            reverseSequence_.push_back(oligo::getBase((~(*bcl)) & 3, true));
            forwardQuality_.push_back(((unsigned char)*bcl) >> 2);
            reverseQuality_.push_back(((unsigned char)*bcl) >> 2);
        }
        else
        {
            forwardSequence_.push_back('n'); // so that it mismatches with 'N' in the reference
            reverseSequence_.push_back('n'); // so that it mismatches with 'N' in the reference
            forwardQuality_.push_back(2);
            reverseQuality_.push_back(2);
        }
    }
    std::reverse(reverseSequence_.begin(), reverseSequence_.end());
    std::reverse(reverseQuality_.begin(), reverseQuality_.end());
}

std::ostream &operator<<(std::ostream &os, const Read &read)
{
    os << "Read:" << std::endl;
    os << "    "; BOOST_FOREACH(const char c, read.getForwardSequence()) {os << c;} os << std::endl;
    os << "    "; BOOST_FOREACH(const char c, read.getForwardQuality()) {os << c;} os << std::endl;
    os << "    "; BOOST_FOREACH(const char c, read.getReverseSequence()) {os << c;} os << std::endl;
    os << "    "; BOOST_FOREACH(const char c, read.getReverseQuality()) {os << c;} os << std::endl;
    return os;
}

} // namespace alignment
} // namespace isaac
