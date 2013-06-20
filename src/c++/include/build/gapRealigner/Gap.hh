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
 ** \file Gap.hh
 **
 ** Gap realigner implementation details.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_GAP_REALIGNER_GAP_HH
#define iSAAC_BUILD_GAP_REALIGNER_GAP_HH

#include "alignment/Cigar.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace build
{
namespace gapRealigner
{
struct Gap
{
    Gap(const reference::ReferencePosition pos, const int length) :
        pos_(pos), length_(length){}
    /**
     * \brief first position of the indel in reference
     */
    reference::ReferencePosition pos_;
    /**
     * \brief positive value indicates deletion from reference. negative - deletion from data.
     *        zero - just an indication of a position in the reference. Intented to be used for lookups.
     */
    int length_;

    unsigned getLength() const {return std::abs(length_);}
    alignment::Cigar::OpCode getOpCode() const { return isInsertion() ? alignment::Cigar::INSERT : alignment::Cigar::DELETE;}
    bool isInsertion() const {return 0 > length_;}
    bool isDeletion() const {return 0 < length_;}

    bool operator ==(const Gap &right) const
    {
        return pos_ == right.pos_ && length_ == right.length_;
    }

    reference::ReferencePosition getBeginPos() const
    {
        return pos_;
    }
    reference::ReferencePosition getEndPos(const bool fatInserstions) const
    {
        return (isDeletion() || fatInserstions) ? pos_ + std::abs(length_) : pos_;
    }

    reference::ReferencePosition getDeletionEndPos() const
    {
        ISAAC_ASSERT_MSG(isDeletion(), "Expected a deletion gap. Got: " << *this);
        return pos_ + std::abs(length_);
    }

    friend std::ostream &operator << (std::ostream &os, const Gap& gap);
};

typedef std::vector<gapRealigner::Gap> Gaps;

inline std::ostream &operator << (std::ostream &os, const Gap& gap)
{
    return os << "Gap(" << gap.pos_ << "," << gap.length_ << ")";
}

struct GapsRange : public std::pair<Gaps::const_iterator, Gaps::const_iterator>
{
    typedef std::pair<Gaps::const_iterator, Gaps::const_iterator> BaseType;
    GapsRange(Gaps::const_iterator f, Gaps::const_iterator s): BaseType(f, s) {}
    GapsRange (const GapsRange &that) : BaseType(that.first, that.second){}
    GapsRange(){}
    bool empty() const {return second == first;}
    unsigned size() const {return second - first;}
};

inline std::ostream &operator <<(std::ostream &os, const GapsRange &gaps)
{
    if (gaps.second == gaps.first)
    {
        return os << "(no gaps)";
    }

    BOOST_FOREACH(const gapRealigner::Gap &gap, std::make_pair(gaps.first, gaps.second))
    {
        os << gap << ",";
    }
    return os;
}

} //namespace gapRealigner

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_GAP_REALIGNER_GAP_HH
