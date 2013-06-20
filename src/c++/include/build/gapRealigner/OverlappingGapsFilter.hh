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
 ** \file OverlappingGapsFilter.hh
 **
 ** Gap realigner implementation details.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_GAP_REALIGNER_OVERLAPPING_GAPS_FILTER_HH
#define iSAAC_BUILD_GAP_REALIGNER_OVERLAPPING_GAPS_FILTER_HH

#include "build/gapRealigner/Gap.hh"
#include "common/BitHacks.hh"

namespace isaac
{
namespace build
{
namespace gapRealigner
{

class OverlappingGapsFilter
{
    static const unsigned MAX_TRACKED_OVERLAPS = 30;
    static const unsigned MAX_TRACKED_DELETIONS = 30;
    const unsigned maxChoice_;
    typedef common::FiniteCapacityVector<unsigned, MAX_TRACKED_OVERLAPS> Overlaps;
    Overlaps overlappingGaps_;
public:
    OverlappingGapsFilter(const gapRealigner::GapsRange &gapsRange)
    : maxChoice_ (gapsRange.size() > MAX_TRACKED_DELETIONS ? 0 : (1 << gapsRange.size()) - 1),
      overlappingGaps_(maxChoice_ ? findOverlaps(gapsRange) : Overlaps())
    {
//        ISAAC_THREAD_CERR << "OverlappingGapsFilter: " << gapsRange << std::endl;
//        ISAAC_THREAD_CERR << "Overlaps: " << std::endl;
//        BOOST_FOREACH(unsigned r, overlappingGaps_)
//        {
//            ISAAC_THREAD_CERR << "Overlap: " << std::hex << r << std::endl;
//        }
    }

    unsigned overlapsCount() const {return overlappingGaps_.size();}
    unsigned overlap(std::size_t n) const {return overlappingGaps_[n];}

    unsigned findOverlaps(const unsigned combination) const
    {
        BOOST_FOREACH(const unsigned overlap, overlappingGaps_)
        {
            const unsigned overlappingGaps = combination & overlap;
            if (overlappingGaps && 1 < countBitsSet(overlappingGaps))
            {
                return overlappingGaps;
            }
        }
        return 0;
    }

    unsigned next(unsigned combination) const
    {
        unsigned increment = 1;
        while (combination < maxChoice_)
        {
            combination += increment;
            const unsigned overlappingGaps = findOverlaps(combination);
            if (!overlappingGaps)
            {
                return combination;
            }
            increment = 1 << lsbSet(overlappingGaps);
        }
        return 0;
    }
private:
    static Overlaps findOverlaps(const gapRealigner::GapsRange &gapsRange);
};

} // namespace gapRealigner
} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_GAP_REALIGNER_OVERLAPPING_GAPS_FILTER_HH
