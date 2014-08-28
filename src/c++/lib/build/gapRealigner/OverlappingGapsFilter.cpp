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
 ** \file OverlappingGapsFilter.cpp
 **
 ** Gap realigner implementation details.
 **
 ** \author Roman Petrovski
 **/

#include "build/gapRealigner/OverlappingGapsFilter.hh"

namespace isaac
{
namespace build
{
namespace gapRealigner
{

typedef std::pair<unsigned, reference::ReferencePosition> GapEnd;
bool orderByEndPosAndIndex(const GapEnd &left, const GapEnd &right)
{
    return left.second < right.second || (left.second == right.second && left.first < right.first);
}

OverlappingGapsFilter::Overlaps OverlappingGapsFilter::findOverlaps(const gapRealigner::GapsRange &gapsRange)
{
    static const unsigned DELETION_END_INDEX_OFFSET = 0;
    static const unsigned DELETION_START_INDEX_OFFSET = 1024;
    static const unsigned INSERTION_INDEX_OFFSET = 2048;
    common::FiniteCapacityVector<GapEnd, MAX_TRACKED_DELETIONS*2> gapEnds;

    int gapIndex = 0;
    for (GapsRange::first_type it = gapsRange.first; gapsRange.second != it; ++it, ++gapIndex)
    {
        const Gap &gap = *it;
        if (gap.isDeletion())
        {
            // for deletions both start end end is stored. Deletion ends are in the lowest range so that
            // they get processed before anything else.
            gapEnds.push_back(GapEnd(gapIndex + DELETION_START_INDEX_OFFSET, gap.getBeginPos()));
            gapEnds.push_back(GapEnd(gapIndex + DELETION_END_INDEX_OFFSET, gap.getEndPos(false)));
        }
        else
        {
            // For insertions only the end is stored.
            gapEnds.push_back(GapEnd(gapIndex + INSERTION_INDEX_OFFSET, gap.getEndPos(false)));
        }
    }
    std::sort(gapEnds.begin(), gapEnds.end(), orderByEndPosAndIndex);

    Overlaps ret;
    // bitmask of insertions that are open. Bits cleared whenever currentPos changes
    unsigned lastInsertionMask = 0;
    reference::ReferencePosition lastInsertionPos;
    unsigned openDeletions = 0;
    unsigned openInsertions = 0;

    ret.push_back(0);
    bool lastWasDeletionClose = true;
    BOOST_FOREACH(const GapEnd &gapEnd, gapEnds)
    {
        if (DELETION_START_INDEX_OFFSET > gapEnd.first)
        {
//            ISAAC_THREAD_CERR << "Deletion close: " << gapEnd.second << std::endl;
            // Close deletion
            const unsigned gapIndex = gapEnd.first;
            const unsigned gapMask = (1 << gapIndex);
            if (lastWasDeletionClose)
            {
                ret.back() &= ~gapMask;
            }
            else if (openDeletions + openInsertions > 1)
            {
                ret.push_back(ret.back()  & ~lastInsertionMask & ~gapMask);
                lastInsertionMask = 0;
                openInsertions = 0;
            }
            else
            {
                ISAAC_ASSERT_MSG(1 == openDeletions, "To close a deletion there must be one open");
                // closing single deletion. Just forget about it.
                ret.back() = 0;
            }
            lastWasDeletionClose = true;
            --openDeletions;
        }
        else if (INSERTION_INDEX_OFFSET > gapEnd.first)
        {
//            ISAAC_THREAD_CERR << "Deletion open: " << gapEnd.second << std::endl;
            // Open deletion
            const unsigned gapIndex = gapEnd.first - DELETION_START_INDEX_OFFSET;
            const unsigned gapMask = (1 << gapIndex);
            if (lastInsertionMask && lastInsertionPos != gapEnd.second)
            {
                if (openDeletions + openInsertions > 1)
                {
                    ret.push_back((ret.back() & ~lastInsertionMask) | gapMask);
                }
                else
                {
                    // one or no open non-overlapping insertions, forget about it
                    ret.back() = gapMask;
                }
                lastInsertionMask = 0;
                openInsertions = 0;
            }
            else
            {
                ret.back() |= gapMask;
            }

            ++openDeletions;
            lastWasDeletionClose = false;
        }
        else
        {
//            ISAAC_THREAD_CERR << "Insertion: " << gapEnd.second << std::endl;
            const unsigned gapIndex = gapEnd.first - INSERTION_INDEX_OFFSET;
            const unsigned gapMask = (1 << gapIndex);
            // Add insertion
            if (lastInsertionMask && lastInsertionPos != gapEnd.second)
            {
                if (openDeletions + openInsertions > 1)
                {
                    ret.push_back((ret.back() & ~lastInsertionMask) | gapMask);
                }
                else
                {
                    // one or no open non-overlapping insertions, forget about it
                    ret.back() = gapMask;
                }
                lastInsertionMask = gapMask;
                openInsertions = 1;
            }
            else
            {
                ret.back() |= gapMask;
                lastInsertionMask |= gapMask;
                ++openInsertions;
            }
            lastInsertionPos = gapEnd.second;
            lastWasDeletionClose = false;
        }
    }

    if (openDeletions + openInsertions <= 1)
    {
        ret.pop_back();
    }

    return ret;
}


} // namespace gapRealigner
} // namespace build
} // namespace isaac
