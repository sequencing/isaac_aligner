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
 ** \file ThreadStats.hh
 **
 ** \brief Collects the match finding statistics for one thread.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_FINDER_THREAD_STATS_HH
#define iSAAC_ALIGNMENT_MATCH_FINDER_THREAD_STATS_HH

//#include <string>
//#include <vector>
//#include <boost/format.hpp>
//#include <boost/filesystem.hpp>
//#include <boost/lambda/lambda.hpp>
//#include <boost/lambda/bind.hpp>
//#include <boost/noncopyable.hpp>
//#include <boost/thread/mutex.hpp>
//
#include "alignment/SeedMetadata.hh"
//#include "alignment/SeedId.hh"
#include "alignment/Seed.hh"
//#include "alignment/MatchTally.hh"
//#include "alignment/MatchDistribution.hh"
//#include "common/Threads.hpp"
#include "flowcell/TileMetadata.hh"
//#include "flowcell/ReadMetadata.hh"
//#include "io/BclReader.hh"
//#include "io/FileBufCache.hh"
//#include "io/MatchWriter.hh"
//#include "oligo/Kmer.hh"
//#include "oligo/Permutations.hh"
//#include "reference/ReferenceKmer.hh"
#include "statistics/MatchFinderTileStats.hh"

namespace isaac
{
namespace alignment
{
namespace matchFinder
{

class ThreadStats : public std::vector<statistics::MatchFinderTileStats>
{
public:
    typedef std::vector<statistics::MatchFinderTileStats> BaseType;
    const unsigned seeds_;

    ThreadStats(
        const std::vector<SeedMetadata> &seedMetadataList,
        const flowcell::TileMetadataList &tileMetadataList) :
            BaseType(seedMetadataList.size() * (tileMetadataList.back().getIndex() + 1)),
            seeds_(seedMetadataList.size())
            {}

    const statistics::MatchFinderTileStats &getSeedTileStat(
        const SeedMetadata& seed,
        const flowcell::TileMetadata& tile) const
    {
        return at(index(seed, tile));
    }

    ThreadStats(unsigned seeds, size_t size) :
        BaseType(size), seeds_(seeds) {}

    void recordMatch(const Seed &seed, const unsigned count, const unsigned repeatThreshold)
    {
        if (0 == count)
        {
            recordNoMatch(seed);
        }
        else if (1 == count)
        {
            recordUniqueMatch(seed);
        }
        else if (count < repeatThreshold)
        {
            recordRepeatMatch(seed, count);
        }
        else
        {
            recordTooManyRepeats(seed);
        }
    }
    void recordUniqueMatch(const Seed &seed) {
        ++at(index(seed)).uniqueMatchSeeds_;
    }

    void recordNoMatch(const Seed &seed) {
        ++at(index(seed)).noMatchSeeds_;
    }

    void recordRepeatMatch(const Seed &seed, const unsigned count) {
        ++at(index(seed)).repeatMatchSeeds_;
        ++at(index(seed)).repeatMatches_ += count;
    }

    void recordTooManyRepeats(const Seed &seed) {
        ++at(index(seed)).tooManyRepeatsSeeds_;
    }

    ThreadStats & operator =(const ThreadStats &that) {
        assert(that.seeds_ == seeds_);
        assert(that.size() == size());
        BaseType &base(*this);
        base = that;
        return *this;
    }

    const ThreadStats operator +(const ThreadStats &right) const
    {
        assert(right.seeds_ == seeds_);
        assert(right.size() == size());
        ThreadStats ret(seeds_, size());
        std::transform(begin(), end(), right.begin(), ret.begin(),
                       std::plus<statistics::MatchFinderTileStats>());
        return ret;
    }

private:
    unsigned index(const Seed &seed) const {
//        std::cerr << "index for seed id: " << seed << " = " << seed.getSeedId().getTile() * seeds_ + seed.getSeedId().getSeed() << "\n";
        return seed.getSeedId().getTile() * seeds_ + seed.getSeedId().getSeed();
    }
    unsigned index(const SeedMetadata& seed, const flowcell::TileMetadata& tile) const {
//        std::cerr << "index for seed metatdata and tile: " << seed << " " << tile << " = " << tile.getIndex() * seeds_ + seed.getIndex() << "\n";
        return tile.getIndex() * seeds_ + seed.getIndex();
    }

    static unsigned getMaxTileIndex(const flowcell::TileMetadataList &tileMetadataList)
    {
        return std::max_element(tileMetadataList.begin(), tileMetadataList.end(),
                                boost::bind(&flowcell::TileMetadata::getIndex, _1)<
                                boost::bind(&flowcell::TileMetadata::getIndex, _2))->getIndex();
//        return std::accumulate(tileMetadataList.begin(), tileMetadataList.end(),
//                               0, boost::bind(&std::max<unsigned>, _1,
//                                              boost::bind(&flowcell::TileMetadata::getIndex, _2)));
    }

};

} // namespace matchFinder
} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_FINDER_THREAD_STATS_HH
