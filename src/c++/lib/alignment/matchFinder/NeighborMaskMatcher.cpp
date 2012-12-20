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
 ** \file NeighborMaskMatcher.cpp
 **
 ** Top Matches with up to 4 mismatches in the seed suffix.
 **
 ** \author Come Raczy
 **/

#include <boost/foreach.hpp>

#include "alignment/matchFinder/NeighborMaskMatcher.hh"

namespace isaac
{
namespace alignment
{
namespace matchFinder
{

inline bool areNeighbors(oligo::Kmer lhs, oligo::Kmer rhs, const unsigned distance)
{
    ISAAC_ASSERT_MSG(lhs != rhs, "Neighbor matcher is really inefficient way to catch exact matches. "
        "It also can miss some. Exact matches must be handled elsewhere.");

    const unsigned suffixBits = oligo::kmerBitLength / 2;
    ISAAC_ASSERT_MSG((lhs >> suffixBits) == (rhs >> suffixBits), "Prefix must match for neighbor matching");

    unsigned shiftedBits = 0;
    unsigned mismatchCount = 0;
    while (distance >= mismatchCount && shiftedBits < suffixBits)
    {
        if ((lhs & 3) != (rhs & 3))
        {
            ++mismatchCount;
        }
        lhs >>= 2;
        rhs >>= 2;
        shiftedBits += 2;
    }
    return (distance >= mismatchCount);
}

/**
 * \brief For every seed of read read writes out a no-mach so that MatchSelector can properly
 *        count the total number of clusters.
 */
void NeighborMaskMatcher::generateNoMatches(
    const SeedIterator currentSeed,
    const SeedIterator nextSeed,
    io::TileMatchWriter &matchWriter)
{
    for (SeedIterator seed = currentSeed; nextSeed > seed; ++seed)
    {
        if (!foundMatches_.isReadComplete(seed->getTile(), seed->getCluster(), seedMetadataList_[seed->getSeedIndex()].getReadIndex()))
        {
            // closing read here is fine as this is the last step of match generation and the only thing checking
            // for whether the reads is closed is the line above.
            foundMatches_.markReadComplete(seed->getTile(), seed->getCluster(), seedMetadataList_[seed->getSeedIndex()].getReadIndex());
            matchWriter.write(*seed, reference::ReferencePosition(reference::ReferencePosition::NoMatch));
        }
    }
}

/**
 * \brief For every seed a read writes out a too-many-mach so that MatchSelector can properly
 *        count the total number of clusters and ignore the seeds that match too much when picking the best template.
 */
void NeighborMaskMatcher::generateTooManyMatches(
    const SeedIterator currentSeed,
    const SeedIterator nextSeed,
    io::TileMatchWriter &matchWriter)
{
    if (ignoreRepeats_)
    {
        generateNoMatches(currentSeed, nextSeed, matchWriter);
    }
    else
    {
        for (SeedIterator seed = currentSeed; nextSeed > seed; ++seed)
        {
            foundMatches_.markReadComplete(seed->getTile(), seed->getCluster(), seedMetadataList_[seed->getSeedIndex()].getReadIndex());
            matchWriter.write(*seed, reference::ReferencePosition(reference::ReferencePosition::TooManyMatch));
        }
    }
}


void NeighborMaskMatcher::matchNeighborsMask(
    const std::vector<Seed>::const_iterator beginSeeds,
    const std::vector<Seed>::const_iterator endSeeds,
    const oligo::Kmer mask,
    const unsigned maskBases,
    MatchDistribution &matchDistribution,
    std::vector<reference::ReferenceKmer> &threadRepeatList,
    std::vector<reference::ReferenceKmer> &threadNeighborsList,
    io::TileMatchWriter &matchWriter,
    std::istream &reference)
{
    const clock_t start = clock();
    ISAAC_THREAD_CERR << "Finding neighbors matches for mask " <<
        oligo::Bases<oligo::BITS_PER_BASE, oligo::Kmer>(mask, maskBases) << std::endl;
    unsigned matchCounters[] = {0, 0};
    // all memory reservation must have been done outside the threaded code
    assert(threadNeighborsList.capacity() >= repeatThreshold_ + 1);
    threadNeighborsList.clear();
    SeedIterator nextSeed = beginSeeds;
    reference::ReferenceKmer nextReference;
    char *readBuffer = reinterpret_cast<char *>(&nextReference);
    reference.read(readBuffer, sizeof(nextReference));
    const unsigned suffixBits = oligo::kmerBitLength / 2;
    oligo::Kmer currentPrefix = 0;
    if (reference)
    {
        while(endSeeds != nextSeed)
        {
            // Load the list of reference positions matching the current prefix, if needed
            if ((beginSeeds == nextSeed) || (currentPrefix != (nextSeed->getKmer() >> suffixBits)))
            {
                if (!reference)
                {
                    // we need a new prefix from the reference, but the refernce is out of data.
                    // this is where we stop then.
                    break;
                }
                threadNeighborsList.clear();
                currentPrefix = nextSeed->getKmer() >> suffixBits;
                // skip the reference position where the prefix is too small
                while (reference && (currentPrefix > (nextReference.getKmer() >> suffixBits)))
                {
                    reference.read(readBuffer, sizeof(nextReference));
                }
                // store the reference positions where the prefix equals currentPrefix
                while (reference && currentPrefix == (nextReference.getKmer() >> suffixBits))
                {
                    // all memory reservation must have been done outside the threaded code
                    assert(threadNeighborsList.capacity() > threadNeighborsList.size());
                    if (threadNeighborsList.size() < neighborhoodSizeThreshold_)
                    {
                        threadNeighborsList.push_back(reference::ReferenceKmer(nextReference.getKmer(),
                                                                               nextReference.getTranslatedPosition(contigKaryotypes_)));
                    }
                    reference.read(readBuffer, sizeof(nextReference));
                }
            }
            // identify all the seeds with the same k-mer
            const SeedIterator currentSeed = nextSeed;
            while ((endSeeds != nextSeed) && (currentSeed->getKmer() == nextSeed->getKmer()))
            {
                ++nextSeed;
            }

            // skip the neighborhood search if there are too many reference positions with this prefix
            if (threadNeighborsList.size() >= neighborhoodSizeThreshold_)
            {
                generateNoMatches(currentSeed, nextSeed, matchWriter);
                continue;
            }

            // Generate the list of reference positions within a distance of 4 to the currentSeed k-mer
            threadRepeatList.clear();
            std::vector<reference::ReferenceKmer>::const_iterator currentNeighbor = threadNeighborsList.begin();
            while ((threadNeighborsList.end() != currentNeighbor) && (threadRepeatList.size() < repeatThreshold_))
            {
                // exact matches have already been handled by exact matcher
                if (currentNeighbor->getKmer() != currentSeed->getKmer())
                {
                    // the more inexact matching is allowed, the more misplaced reads we're getting
                    if (areNeighbors(currentNeighbor->getKmer(), currentSeed->getKmer(), 1))
                    {
                        if (currentNeighbor->getReferencePosition().isTooManyMatch())
                        {
                            // for high repeat hits ensure TooManyMatches gets generated down below
                            threadRepeatList.resize(repeatThreshold_);
                            break;
                        }
                        else
                        {
                            threadRepeatList.push_back(*currentNeighbor);
                        }
                    }
                }
                ++currentNeighbor;
            }
            if(threadRepeatList.empty())
            {
                generateNoMatches(currentSeed, nextSeed, matchWriter);
            }
            else if (threadRepeatList.size() >= repeatThreshold_)
            {
                // MatchSelector is smart enough to filter out random neighbor matches.
                // Silencing good seed matches this way actually causes trouble.
                generateTooManyMatches(currentSeed, nextSeed, matchWriter);
            }
            else
            {
                // generate the matches for each seed if needed
                for (SeedIterator seed = currentSeed; nextSeed > seed; ++seed)
                {
                    matchCounters[seed->isReverse()] += threadRepeatList.size();
                    BOOST_FOREACH(const reference::ReferenceKmer &repeat, threadRepeatList)
                    {
                        foundMatches_.markReadComplete(seed->getTile(), seed->getCluster(), seedMetadataList_[seed->getSeedIndex()].getReadIndex());
                        reference::ReferencePosition repeatNeighborMatchPosition = repeat.getReferencePosition();
                        // Because neighbor mask matcher finds only half of the inexact matches, cases
                        // where reference reverse-compliments itself can be particularly misleading:
                        // It will find one location but not the other.
                        // force neighbor bit so that MatchSelector gets bit more cautious about these matches.
                        repeatNeighborMatchPosition.setNeighbors(true);
                        matchWriter.write(*seed, repeatNeighborMatchPosition);
                    }
                }
                // assume that the repeat resolution will be reasonably uniform
                const size_t matchesPerRepeat = size_t(round(static_cast<double>(nextSeed - currentSeed) / static_cast<double>(threadRepeatList.size())));
                BOOST_FOREACH(const reference::ReferenceKmer &repeat, threadRepeatList)
                {
                    const reference::ReferencePosition position = repeat.getReferencePosition();
                    matchDistribution.addMatches(position.getContigId(), position.getPosition(), matchesPerRepeat);
                }
            }

        }
    }
    ISAAC_THREAD_CERR << "Finding neighbors matches done in " << (clock() - start) / 1000 << " ms for mask " <<
        oligo::Bases<oligo::BITS_PER_BASE, oligo::Kmer>(mask, maskBases) << std::endl;

    ISAAC_THREAD_CERR <<
        "Found " << (matchCounters[0] + matchCounters[1]) <<
        " matches (" << matchCounters[0] << " forward, " << matchCounters[1] << " reverse)"
        " for mask " << oligo::Bases<oligo::BITS_PER_BASE, oligo::Kmer>(mask, maskBases) <<
        " and " << (endSeeds - beginSeeds) << " kmers"
        " in range [" << oligo::Bases<oligo::BITS_PER_BASE, oligo::Kmer>(beginSeeds->getKmer(), oligo::kmerLength) <<
        "," << std::setbase(16) << (beginSeeds == endSeeds ?
            oligo::Bases<oligo::BITS_PER_BASE, oligo::Kmer>(beginSeeds->getKmer(), oligo::kmerLength) :
            oligo::Bases<oligo::BITS_PER_BASE, oligo::Kmer>((endSeeds - 1)->getKmer(), oligo::kmerLength)) <<
        "]" << std::endl;
}

} // namespace matchFinder
} // namespace alignment
} // namespace isaac
