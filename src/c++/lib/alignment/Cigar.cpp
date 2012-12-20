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
 ** \file Cigar.cpp
 **
 ** \brief See Cigar.hh
 ** 
 ** \author Come Raczy
 **/

#include "alignment/Cigar.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

std::string Cigar::toString() const
{
    return toString(begin(), end());
}

std::string Cigar::toString(const unsigned offset, const unsigned length) const
{
    ISAAC_ASSERT_MSG(this->size() >= offset + length, "Requested end is outside of cigarBuffer");
    return toString(begin() + offset, begin() + offset + length);
}

std::string Cigar::toString(const std::vector<uint32_t> &cigarBuffer, unsigned offset, unsigned length)
{
    ISAAC_ASSERT_MSG(cigarBuffer.size() >= offset + length, "Requested end is outside of cigarBuffer");
    return toString(cigarBuffer.begin() + offset, cigarBuffer.begin() + offset + length);
}

} // namespace alignment
} // namespace isaac
