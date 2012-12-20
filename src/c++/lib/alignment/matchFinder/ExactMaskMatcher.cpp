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
 ** \file ExactMaskMatcher.cpp
 **
 ** Top level component to perform an alignment.
 **
 ** \author Come Raczy
 **/

#include <boost/foreach.hpp>

#include "alignment/matchFinder/ExactMaskMatcher.hh"

namespace isaac
{
namespace alignment
{
namespace matchFinder
{

///// find the next range of non N-masked seeds after nextSeed
//void advance(SeedIterator &currentSeed, SeedIterator &nextSeed, const SeedIterator endSeeds)
//{
//    currentSeed = nextSeed;
//    while ((endSeeds != currentSeed) && currentSeed->getNMask())
//    {
//        ++currentSeed;
//    }
//    nextSeed = currentSeed;
//    while ((endSeeds != nextSeed) && (currentSeed->getKmer() == nextSeed->getKmer()))
//    {
//        ++nextSeed;
//    }
//}
/**
 * \brief For every seed of read read writes out a no-mach so that MatchSelector can properly
 *        count the total number of clusters.
 */
void ExactMaskMatcher::generateNoMatches(
    const SeedIterator currentSeed,
    const SeedIterator nextSeed,
    io::TileMatchWriter &matchWriter)
{
    for (SeedIterator seed = currentSeed; nextSeed > seed; ++seed)
    {
        if (!foundExactMatchesOnly_.isReadComplete(seed->getTile(), seed->getCluster(), seedMetadataList_[seed->getSeedIndex()].getReadIndex()))
        {
            matchWriter.write(*seed, reference::ReferencePosition(reference::ReferencePosition::NoMatch));
        }
    }
}

/**
 * \brief For every seed, writes out a too-many-mach so that MatchSelector can properly
 *        count the total number of clusters and ignore those seeds when picking the best template.
 *        Completes the read so that it does not get misplaced by Neighbor matches
 */
void ExactMaskMatcher::generateTooManyMatches(
    const SeedIterator currentSeed,
    const SeedIterator nextSeed,
    io::TileMatchWriter &matchWriter)
{
    for (SeedIterator seed = currentSeed; nextSeed > seed; ++seed)
    {
        const unsigned readIndex = seedMetadataList_[seed->getSeedIndex()].getReadIndex();

        matchWriter.write(*seed, reference::ReferencePosition(reference::ReferencePosition::TooManyMatch));
        if (closeRepeats_)
        {
            foundExactMatchesOnly_.markReadComplete(seed->getTile(), seed->getCluster(), readIndex);
        }
    }
}

// Initial implementation that works only for exact matches
void ExactMaskMatcher::matchMask(
    const std::vector<Seed>::const_iterator beginSeeds,
    const std::vector<Seed>::const_iterator endSeeds,
    const oligo::Kmer mask,
    const unsigned maskBases,
    MatchDistribution &matchDistribution,
    std::vector<reference::ReferenceKmer> &threadRepeatList,
    io::TileMatchWriter &matchWriter,
    std::istream &reference)
{
    const clock_t start = clock();
    ISAAC_THREAD_CERR << "Finding exact matches for mask " << mask << std::endl;
    unsigned matchCounters[] = {0, 0};
    unsigned repeatCounters[] = {0, 0};
    unsigned highRepeatCounters[] = {0, 0};
    // all memory reservation must have been done outside the threaded code
    assert(threadRepeatList.capacity() >= repeatThreshold_ + 1);
    typedef std::vector<Seed>::const_iterator SeedIterator;
    SeedIterator nextSeed = beginSeeds;
    reference::ReferenceKmer nextReference;
    char *readBuffer = reinterpret_cast<char *>(&nextReference);
    reference.read(readBuffer, sizeof(nextReference));
    while(endSeeds != nextSeed)
    {
        // identify all the seeds with the same k-mer
        const SeedIterator currentSeed = nextSeed;
        while ((endSeeds != nextSeed) && (currentSeed->getKmer() == nextSeed->getKmer()))
        {
            ++nextSeed;
        }
        // discard reference positions with smaller k-mer
        while (reference && currentSeed->getKmer() > nextReference.getKmer())
        {
            reference.read(readBuffer, sizeof(nextReference));
        }
        // Generate the list of reference positions matching the currentSeed
        threadRepeatList.clear();
        while (reference && (currentSeed->getKmer() == nextReference.getKmer()))
        {
            if (threadRepeatList.size() < repeatThreshold_)
            {
                threadRepeatList.push_back(reference::ReferenceKmer(nextReference.getKmer(),
                                                                    nextReference.getTranslatedPosition(contigKaryotypes_)));
            }
            reference.read(readBuffer, sizeof(nextReference));
        }
        // generate the matches for each seed
        if (threadRepeatList.empty())
        {
            if (storeNomatches_)
            {
                generateNoMatches(currentSeed, nextSeed, matchWriter);
            }
        }
        else if (threadRepeatList.size() >= repeatThreshold_)
        {
            generateTooManyMatches(currentSeed, nextSeed, matchWriter);
            for (SeedIterator seed = currentSeed; nextSeed > seed; ++seed)
            {
                repeatCounters[seed->isReverse()] += threadRepeatList.size();
            }
        }
        else
        {
            const reference::ReferencePosition anyPosition = threadRepeatList.front().getReferencePosition();
            if (anyPosition.isTooManyMatch())
            {
                ISAAC_ASSERT_MSG(1 == threadRepeatList.size(), "Reference must contain no more than 1 TooManyMatch entries per kmer");
                generateTooManyMatches(currentSeed, nextSeed, matchWriter);
                for (SeedIterator seed = currentSeed; nextSeed > seed; ++seed)
                {
                    ++highRepeatCounters[seed->isReverse()];
                }
            }
            else
            {
                for (SeedIterator seed = currentSeed; nextSeed > seed; ++seed)
                {
                    const unsigned readIndex = seedMetadataList_[seed->getSeedIndex()].getReadIndex();
                    matchCounters[seed->isReverse()] += threadRepeatList.size();

                    BOOST_FOREACH(const reference::ReferenceKmer &repeat, threadRepeatList)
                    {
                        matchWriter.write(*seed, repeat.getReferencePosition());
                    }
                    if (ignoreNeighbors_ || !anyPosition.hasNeighbors())
                    {
                        foundExactMatchesOnly_.markReadComplete(seed->getTile(), seed->getCluster(), readIndex);
                    }
                }
                // update the MatchDistribution
                // assume that the repeat resolution is reasonably uniform
                const size_t matchesPerRepeat = size_t(round(static_cast<double>(nextSeed - currentSeed) / static_cast<double>(threadRepeatList.size())));
                BOOST_FOREACH(const reference::ReferenceKmer &repeat, threadRepeatList)
                {
                    const reference::ReferencePosition position = repeat.getReferencePosition();
                    matchDistribution.addMatches(position.getContigId(), position.getPosition(), matchesPerRepeat);
                }
            }
        }
    }

    ISAAC_THREAD_CERR << "Finding exact matches done in " << (clock() - start) / 1000 << " ms for mask " <<
        mask << std::endl;

    ISAAC_THREAD_CERR <<
        "Found " << (matchCounters[0] + matchCounters[1]) <<
        " matches (" << matchCounters[0] << " forward, " << matchCounters[1] << " reverse),"
        " repeats at least (" << repeatCounters[0] << " forward, " << repeatCounters[1] << " reverse),"
        " high repeats (" << highRepeatCounters[0] << " forward, " << highRepeatCounters[1] << " reverse)"
        " for mask " << mask <<
        " and " << (endSeeds - beginSeeds) << " kmers"
        " in range [" << oligo::Bases<oligo::BITS_PER_BASE, oligo::Kmer>(beginSeeds->getKmer(), oligo::kmerLength) <<
        "," << (beginSeeds == endSeeds ?
            oligo::Bases<oligo::BITS_PER_BASE, oligo::Kmer>(beginSeeds->getKmer(), oligo::kmerLength) :
            oligo::Bases<oligo::BITS_PER_BASE, oligo::Kmer>((endSeeds - 1)->getKmer(), oligo::kmerLength)) <<
        "]" << std::endl;

}

} // namespace matchFinder
} // namespace alignment
} // namespace isaac
