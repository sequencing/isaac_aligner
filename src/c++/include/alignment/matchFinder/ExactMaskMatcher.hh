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

template <typename KmerT>
class ExactMaskMatcher: boost::noncopyable
{
    typedef Seed<KmerT> SeedT;
    typedef typename std::vector<SeedT>::const_iterator SeedIterator;
    typedef reference::ReferenceKmer<KmerT> ReferenceKmerT;

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
        const SeedIterator beginSeeds,
        const SeedIterator endSeeds,
        const unsigned mask,
        MatchDistribution &matchDistribution,
        std::vector<ReferenceKmerT> &threadRepeatList,
        io::TileMatchWriter &matchWriter,
        std::istream &reference);

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
