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
 ** \file ExactMaskMatcher.hh
 **
 ** \brief Walks along the sorted seeds and sorted reference and produces the exact kmer matches.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_FINDER_NEIGHBOR_MASK_MATCHER_HH
#define iSAAC_ALIGNMENT_MATCH_FINDER_NEIGHBOR_MASK_MATCHER_HH

#include <vector>

#include "alignment/MatchDistribution.hh"
#include "alignment/matchFinder/ThreadStats.hh"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "alignment/Seed.hh"
#include "io/MatchWriter.hh"
#include "oligo/Kmer.hh"
#include "reference/ReferenceKmer.hh"

namespace isaac
{
namespace alignment
{
namespace matchFinder
{

class NeighborMaskMatcher: boost::noncopyable
{
public:
    NeighborMaskMatcher(
        const bool ignoreRepeats,
        const unsigned repeatThreshold,
        const unsigned neighborhoodSizeThreshold,
        const SeedMetadataList &seedMetadataList,
        const std::vector<unsigned>& contigKaryotypes,
        matchFinder::TileClusterInfo &foundExactMatchesOnly) :
        ignoreRepeats_(ignoreRepeats), repeatThreshold_(repeatThreshold), neighborhoodSizeThreshold_(neighborhoodSizeThreshold),
        seedMetadataList_(seedMetadataList), contigKaryotypes_(contigKaryotypes),
        foundMatches_(foundExactMatchesOnly){}
    /// walks along the sorted seeds and sorted reference and produces the matches
    void matchNeighborsMask(
        const std::vector<Seed>::const_iterator beginSeeds,
        const std::vector<Seed>::const_iterator endSeeds,
        const oligo::Kmer mask,
        const unsigned maskBases,
        MatchDistribution &matchDistribution,
        std::vector<reference::ReferenceKmer> &threadRepeatList,
        std::vector<reference::ReferenceKmer> &threadNeighborsList,
        io::TileMatchWriter &matchWriter,
        std::istream &reference);

private:
    const bool ignoreRepeats_;
    const unsigned repeatThreshold_;
    const unsigned neighborhoodSizeThreshold_;
    const SeedMetadataList &seedMetadataList_;
    const std::vector<unsigned>& contigKaryotypes_;

    matchFinder::TileClusterInfo &foundMatches_;

    typedef std::vector<Seed>::const_iterator SeedIterator;
    void generateNoMatches(
        const SeedIterator currentSeed,
        const SeedIterator nextSeed,
        io::TileMatchWriter &matchWriter);
    void generateTooManyMatches(
        const SeedIterator currentSeed,
        const SeedIterator nextSeed,
        io::TileMatchWriter &matchWriter);
};

} // namespace matchFinder
} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_FINDER_NEIGHBOR_MASK_MATCHER_HH
