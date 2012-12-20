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

#ifndef iSAAC_ALIGNMENT_MATCH_FINDER_EXACT_MASK_MATCHER_HH
#define iSAAC_ALIGNMENT_MATCH_FINDER_EXACT_MASK_MATCHER_HH

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

class ExactMaskMatcher: boost::noncopyable
{
    const bool closeRepeats_;
    const bool storeNomatches_;
    const unsigned repeatThreshold_;
    const bool ignoreNeighbors_;
    const std::vector<unsigned>& contigKaryotypes_;
    const SeedMetadataList &seedMetadataList_;

    matchFinder::TileClusterInfo &foundExactMatchesOnly_;

public:
    ExactMaskMatcher(
        const bool closeRepeats,
        const bool storeNomatches,
        const unsigned repeatThreshold,
        const bool ignoreNeighbors,
        const SeedMetadataList &seedMetadataList,
        const std::vector<unsigned>& contigKaryotypes,
        matchFinder::TileClusterInfo &foundExactMatchesOnly) :
        closeRepeats_(closeRepeats), storeNomatches_(storeNomatches), repeatThreshold_(repeatThreshold),
        ignoreNeighbors_(ignoreNeighbors), contigKaryotypes_(contigKaryotypes), seedMetadataList_(seedMetadataList),
        foundExactMatchesOnly_(foundExactMatchesOnly){}
    /// walks along the sorted seeds and sorted reference and produces the matches
    void matchMask(
        const std::vector<Seed>::const_iterator beginSeeds,
        const std::vector<Seed>::const_iterator endSeeds,
        const oligo::Kmer mask,
        const unsigned maskBases,
        MatchDistribution &matchDistribution,
        std::vector<reference::ReferenceKmer> &threadRepeatList,
        io::TileMatchWriter &matchWriter,
        std::istream &reference);

    typedef std::vector<Seed>::const_iterator SeedIterator;
    void generateTooManyMatches(
        const SeedIterator currentSeed,
        const SeedIterator nextSeed,
        io::TileMatchWriter &matchWriter);
    void generateNoMatches(
        const SeedIterator currentSeed,
        const SeedIterator nextSeed,
        io::TileMatchWriter &matchWriter);
};

} // namespace matchFinder
} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_FINDER_EXACT_MASK_MATCHER_HH
