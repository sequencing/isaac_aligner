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
 ** \file MatchFinderTileStats.hh
 **
 ** \brief MatchFinder statistics helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_STATISTICS_MATCH_FINDER_TILE_STATS_H
#define ISAAC_STATISTICS_MATCH_FINDER_TILE_STATS_H

namespace isaac
{
namespace statistics{

struct MatchFinderTileStats
{
    MatchFinderTileStats() :
        uniqueMatchSeeds_(0), noMatchSeeds_(0), repeatMatchSeeds_(0), tooManyRepeatsSeeds_(0), repeatMatches_(0){}
    unsigned long uniqueMatchSeeds_;
    unsigned long noMatchSeeds_;
    unsigned long repeatMatchSeeds_;
    unsigned long tooManyRepeatsSeeds_;
    unsigned long repeatMatches_;
};

inline const MatchFinderTileStats operator +(MatchFinderTileStats left, const MatchFinderTileStats &right)
{
    left.uniqueMatchSeeds_ += right.uniqueMatchSeeds_;
    left.noMatchSeeds_ += right.noMatchSeeds_;
    left.repeatMatchSeeds_ += right.repeatMatchSeeds_;
    left.tooManyRepeatsSeeds_ += right.tooManyRepeatsSeeds_;
    left.repeatMatches_ += right.repeatMatches_;
    return left;
}

} //namespace statistics
} //namespace isaac

#endif //ISAAC_STATISTICS_MATCH_FINDER_TILE_STATS_H
