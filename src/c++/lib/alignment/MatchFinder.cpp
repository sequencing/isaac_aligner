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
 ** \file MatchFinder.cpp
 **
 ** Top level component to perform an alignment.
 **
 ** \author Come Raczy
 **/

#include <fcntl.h>

#include <algorithm>
#include <numeric>
#include <cassert>
#include <cstring>
#include <cerrno>
#include <ctime>

#include <boost/filesystem.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/system/error_code.hpp>
#include <boost/thread.hpp>

#include "alignment/MatchFinder.hh"
#include "alignment/matchFinder/ExactMaskMatcher.hh"
#include "alignment/matchFinder/NeighborMaskMatcher.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/SystemCompatibility.hh"
#include "common/Threads.hpp"
#include "flowcell/Layout.hh"
#include "io/MatchWriter.hh"
#include "statistics/MatchFinderStatsXml.hh"

namespace isaac
{
namespace alignment
{

static const std::string ABCD("ABCD");

std::vector<std::vector<unsigned> > getReferenceContigKaryotypes(
    const reference::SortedReferenceXmlList &sortedReferenceList)
{
    std::vector<std::vector<unsigned> > ret(sortedReferenceList.size());
    std::vector<std::vector<unsigned> >::iterator it = ret.begin();
    BOOST_FOREACH(const reference::SortedReferenceXml &sortedReference, sortedReferenceList)
    {
        const reference::SortedReferenceXml::Contigs &contigs = sortedReference.getContigs();
        it->reserve(contigs.size());
        std::transform(contigs.begin(), contigs.end(), std::back_inserter(*it),
                       boost::bind(&reference::SortedReferenceXml::Contig::karyotypeIndex_, _1));

        ++it;
    }
    return ret;
}

/**
 * \param   maxRepeatThreshold the maximum number of repeats that will be ever used with the match call.
 *          Required for upfront memory reservation.
 */
MatchFinder::MatchFinder(
    const reference::SortedReferenceXmlList &sortedReferenceList,
    const bfs::path & tempDirectory,
    const flowcell::TileMetadataList &tiles,
    const flowcell::ReadMetadataList &readMetadataList,
    const SeedMetadataList &seedMetadataList,
    const unsigned iteration,
    const bool ignoreNeighbors,
    const bool ignoreRepeats,
    const unsigned repeatThreshold,
    const unsigned neighborhoodSizeThreshold,
    MatchTally &matchTally,
    matchFinder::TileClusterInfo &foundMatches,
    common::ThreadVector &threads,
    const unsigned coresMax,
    const unsigned tempSaversMax,
    const unsigned unavailableFileHandles)
    : kmerSourceMetadataList_(getMaskFilesList(sortedReferenceList, ABCD))
    , referenceContigKaryotypes_(getReferenceContigKaryotypes(sortedReferenceList))
    , seedMetadataList_(seedMetadataList)
    , iteration_(iteration)
    , ignoreNeighbors_(ignoreNeighbors)
    , ignoreRepeats_(ignoreRepeats)
    , repeatThreshold_(repeatThreshold)
    , neighborhoodSizeThreshold_(neighborhoodSizeThreshold)
    , foundExactMatchesOnly_(foundMatches)
    , threads_(threads)
    , threadsMax_(coresMax)
    , maxTilesAtATime_(verifyMaxTileCount(unavailableFileHandles, tempSaversMax, tiles.size()))
    , threadRepeatLists_(threadsMax_, std::vector<reference::ReferenceKmer>(repeatThreshold_ + 1))
    , threadNeighborsLists_(threadsMax_, std::vector<reference::ReferenceKmer>(neighborhoodSizeThreshold_ + 1))
    , threadMatchDistributions_(threadsMax_, MatchDistribution(sortedReferenceList))
    , matchWriter_(matchTally, maxTilesAtATime_, tiles.back().getIndex())
    , threadReferenceFileBuffers_(
        threadsMax_,
        io::FileBufCache<io::FileBufWithReopen>(1, std::ios_base::binary|std::ios_base::in,
                                                std::max_element(kmerSourceMetadataList_.begin(), kmerSourceMetadataList_.end(),
                                                                 boost::bind(&KmerSourceMetadata::getPathSize, _1)<
                                                                 boost::bind(&KmerSourceMetadata::getPathSize, _2))->getPathSize()))
{
    ISAAC_THREAD_CERR << "Constructing the match finder" << std::endl;

    ISAAC_THREAD_CERR << "Constructing the match finder done" << std::endl;
}

/*
 * \brief We have to keep an open output file for each tile in process.
 *        Limit the number of tiles so that we don't go over the ulimit -n
 */
unsigned MatchFinder::verifyMaxTileCount(
    const unsigned unavailableFileHandlesCount,
    const unsigned maxSavers,
    const unsigned tilesCount) const
{
    const unsigned requiredInputFilesPerThread(1); // reference kmers
    const unsigned requiredOutputFilesPerTile(1);  // this implementation stores all matches of a tile in a single file
                                                   // the thread synchronization is used to prevent file corruption

    const unsigned maxTileCount =
        (common::getMaxOpenFiles() -
            unavailableFileHandlesCount -
            requiredInputFilesPerThread * threadsMax_) / requiredOutputFilesPerTile;

    if (0 == maxTileCount) {
        BOOST_THROW_EXCEPTION(common::PreConditionException(
            (boost::format("MatchFinder needs to open at least %d files to process this data set on %d threads."
                " Reduce thread count or increase file handle limit (current is %d)") %
                (threadsMax_ * requiredOutputFilesPerTile) % common::getMaxOpenFiles()).str()));
    }

    const unsigned ret = std::min<unsigned>(tilesCount, std::min(maxSavers, maxTileCount));

    ISAAC_THREAD_CERR << "verifyMaxTileCount returning " << maxTileCount << std::endl;

    return ret;
}

const std::vector<MatchFinder::KmerSourceMetadata> MatchFinder::getMaskFilesList(
    const reference::SortedReferenceXmlList &sortedReferenceList,
    const std::string &permutationName) const
{
    std::vector<KmerSourceMetadata> ret;
    BOOST_FOREACH(const reference::SortedReferenceXml &sortedReference, sortedReferenceList)
    {
        const unsigned maskWidth = sortedReference.getDefaultMaskWidth();
        const unsigned maskCount = oligo::getMaskCount(maskWidth);
        for (oligo::Kmer mask = 0; mask < maskCount; ++mask)
        {
            ret.push_back(
                KmerSourceMetadata(&sortedReference - &sortedReferenceList.front(),
                                   maskWidth, mask,
                                   sortedReference.getMaskFile(permutationName, maskWidth, mask)));
        }
    }
    return ret;
}

void MatchFinder::setTiles(const flowcell::TileMetadataList &tiles)
{
    matchWriter_.reopen(iteration_, tiles);
}

/**
 * \brief Find matches.
 *
 * \param seeds Expects sorted ABCD permutation of seeds. Ns not included.
 */
const std::vector<MatchDistribution> &MatchFinder::findMatches(
    const std::vector<Seed>::iterator seedsBegin,
    const std::vector<std::vector<Seed>::iterator> &referenceSeedBounds,
    const bool findNeighbors,
    const bool finalPass)
{
    return  match(seedsBegin, referenceSeedBounds, findNeighbors, finalPass, ABCD);

    /*
     * TODO: support up to two mismatches. Uncomment the lines below to iterate
     * over all the possible permutations and add the functionality in the
     * matchMask method (see the version that has been commented out) (priority
     * 2).
     */
    //BOOST_FOREACH(const oligo::NamedPermutation &namedPermutation, oligo::permutations)
    //{
    //    oligo::permuteBlocks(namedPermutation.first, seeds);
    //    match(seeds, matchWriter, namedPermutation.second, false);
    //}
    //getBinSizeDistribution(seeds);
}

const std::vector<MatchDistribution> & MatchFinder::match(
    std::vector<Seed>::const_iterator seedsBegin,
    const std::vector<std::vector<Seed>::iterator> &referenceSeedBounds,
    const bool findNeighbors,
    const bool finalPass,
    const std::string &permutationName)
{

    // parallelize the matching across all the masks of all the references.
    KmerSourceMetadataList::const_iterator kmerSourceIterator = kmerSourceMetadataList_.begin();
    threads_.execute(boost::bind(&MatchFinder::matchMaskParallel, this,
                                 boost::ref(seedsBegin),
                                 boost::cref(referenceSeedBounds),
                                 findNeighbors,
                                 finalPass,
                                 boost::ref(kmerSourceIterator),
                                 boost::ref(permutationName),
                                 _1),
                     threadsMax_);

    return threadMatchDistributions_;
}

std::pair<std::vector<Seed>::const_iterator, std::vector<Seed>::const_iterator>
MatchFinder::skipToTheNextMask(
    const std::vector<Seed>::const_iterator currentBegin,
    const std::vector<Seed>::const_iterator seedsEnd,
    const oligo::Kmer currentMask,
    const unsigned maskWidth,
    const bool storeNSeedNoMatches)
{
    const oligo::Kmer endSeed =
        (currentMask << (oligo::kmerBitLength - maskWidth)) | ((~0UL << (oligo::kmerBitLength - maskWidth)) ^ ~0UL);
    // the first possibly-N-seed of the current mask (only the 111.. masks will actually have N-seeds, but it does not matter)
    const Seed firstPossiblyNSeed(endSeed, SMALLEST_N_SEED_ID);
    std::vector<Seed>::const_iterator currentEnd =
        std::lower_bound(currentBegin, seedsEnd, firstPossiblyNSeed, orderByKmerSeedIndex);
    // roll to the first non-N seed so that next thread does not have to deal with our Ns
    std::vector<Seed>::const_iterator nextBegin = std::find_if(currentEnd, seedsEnd, boost::bind(&Seed::isNSeed, _1) != true);

    if (nextBegin != currentEnd)
    {
        ISAAC_THREAD_CERR << "Skipped " << nextBegin - currentEnd << " N-seeds for mask " << currentMask << std::endl;
        if (storeNSeedNoMatches)
        {
            // generate no-match entries for N-containing seeds or else the match selector statistics will never see those clusters
            BOOST_FOREACH(const Seed &seed, std::make_pair(currentEnd, nextBegin))
            {
                // N-seeds don't have a valid seed index. seedMetadataList_[seed.getSeedIndex()] is invalid
//                if (!foundExactMatchesOnly_.isReadComplete(seed.getTile(), seed.getCluster(),
//                                                           seedMetadataList_[seed.getSeedIndex()].getReadIndex()))
                {
                    matchWriter_.write(seed.getSeedId(), reference::ReferencePosition(reference::ReferencePosition::NoMatch));
                }
            }
        }
    }
    return std::make_pair(currentEnd, nextBegin);
}

void MatchFinder::matchMaskParallel(
    std::vector<Seed>::const_iterator &seedsBegin,
    const std::vector<std::vector<Seed>::iterator> &referenceSeedBounds,
    const bool findNeighbors,
    const bool finalPass,
    KmerSourceMetadataList::const_iterator &kmerSourceIterator,
    const std::string &permutationName,
    const unsigned threadNumber)
{
    boost::lock_guard<boost::mutex> lock(mutex_);
    while (kmerSourceIterator < kmerSourceMetadataList_.end() && referenceSeedBounds.back() != seedsBegin)
    {
        KmerSourceMetadataList::const_iterator ourKmerSource = kmerSourceIterator++;
//        const unsigned int maskCount = oligo::getMaskCount(ourKmerSource->maskWidth_);
        const oligo::Kmer currentMask = ourKmerSource->mask_;

        const std::vector<Seed>::const_iterator ourBegin(seedsBegin);
        const std::vector<Seed>::const_iterator seedsEnd(referenceSeedBounds.at(ourKmerSource->referenceIndex_));

        // on the final pass make sure the n-seeds of the open reads get their
        // nomatches stored. Else Match selector stats will report incorrect total cluster count
        std::pair<std::vector<Seed>::const_iterator, std::vector<Seed>::const_iterator> ourEndNextBegin =
            skipToTheNextMask(ourBegin, seedsEnd, currentMask, ourKmerSource->maskWidth_, finalPass);
        seedsBegin = ourEndNextBegin.second;

        const boost::filesystem::path &sortedReferencePath = ourKmerSource->maskFilePath_;
        std::istream threadReferenceFile(threadReferenceFileBuffers_.at(threadNumber).get(sortedReferencePath, io::FileBufWithReopen::SequentialOften));

        {
            common::unlock_guard<boost::mutex> unlock(mutex_);

            if (findNeighbors)
            {
                matchFinder::NeighborMaskMatcher(
                    ignoreRepeats_,
                    repeatThreshold_,  neighborhoodSizeThreshold_,
                    seedMetadataList_,
                    referenceContigKaryotypes_.at(ourKmerSource->referenceIndex_),
                    foundExactMatchesOnly_).matchNeighborsMask(
                        ourBegin, ourEndNextBegin.first, currentMask,
                        ourKmerSource->maskWidth_ / oligo::BITS_PER_BASE,
                        threadMatchDistributions_[threadNumber],
                        threadRepeatLists_[threadNumber],
                        threadNeighborsLists_[threadNumber],
                        matchWriter_,
                        threadReferenceFile);
            }
            else
            {
                matchFinder::ExactMaskMatcher(
                    1 == iteration_, finalPass,
                    repeatThreshold_,
                    ignoreNeighbors_, seedMetadataList_,
                    referenceContigKaryotypes_.at(ourKmerSource->referenceIndex_),
                    foundExactMatchesOnly_).matchMask(
                        ourBegin, ourEndNextBegin.first, currentMask,
                        ourKmerSource->maskWidth_ / oligo::BITS_PER_BASE,
                        threadMatchDistributions_[threadNumber],
                        threadRepeatLists_[threadNumber],
                        matchWriter_,
                        threadReferenceFile);
            }
            if(!threadReferenceFile && !threadReferenceFile.eof())
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to read reference from " + sortedReferencePath.string()));
            }
        }
    }
}

} // namespace alignment
} // namespace isaac
