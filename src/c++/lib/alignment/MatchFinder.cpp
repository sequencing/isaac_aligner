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

std::vector<std::vector<unsigned> > getReferenceContigKaryotypes(
    const reference::SortedReferenceMetadataList &sortedReferenceList)
{
    std::vector<std::vector<unsigned> > ret(sortedReferenceList.size());
    std::vector<std::vector<unsigned> >::iterator it = ret.begin();
    BOOST_FOREACH(const reference::SortedReferenceMetadata &sortedReference, sortedReferenceList)
    {
        const reference::SortedReferenceMetadata::Contigs &contigs = sortedReference.getContigs();
        it->reserve(contigs.size());
        std::transform(contigs.begin(), contigs.end(), std::back_inserter(*it),
                       boost::bind(&reference::SortedReferenceMetadata::Contig::karyotypeIndex_, _1));

        ++it;
    }
    return ret;
}

/**
 * \param   maxRepeatThreshold the maximum number of repeats that will be ever used with the match call.
 *          Required for upfront memory reservation.
 */
template<typename KmerT>
MatchFinder<KmerT>::MatchFinder(
    const reference::SortedReferenceMetadataList &sortedReferenceList,
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
    : kmerSourceMetadataList_(getMaskFilesList(sortedReferenceList))
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
    , threadRepeatLists_(threadsMax_, std::vector<ReferenceKmer>(repeatThreshold_ + 1))
    , threadNeighborsLists_(threadsMax_, std::vector<ReferenceKmer>(neighborhoodSizeThreshold_ + 1))
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
template<typename KmerT>
unsigned MatchFinder<KmerT>::verifyMaxTileCount(
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

template <typename KmerT>
const std::vector<typename MatchFinder<KmerT>::KmerSourceMetadata> MatchFinder<KmerT>::getMaskFilesList(
    const reference::SortedReferenceMetadataList &sortedReferenceList) const
{
    std::vector<KmerSourceMetadata> ret;
    BOOST_FOREACH(const reference::SortedReferenceMetadata &sortedReference, sortedReferenceList)
    {
        BOOST_FOREACH(const reference::SortedReferenceMetadata::MaskFile &mask,
                      sortedReference.getMaskFileList(oligo::KmerTraits<KmerT>::KMER_BASES))
        {
            ret.push_back(KmerSourceMetadata(&sortedReference - &sortedReferenceList.front(),
                                             mask.maskWidth, mask.mask_, mask.path));
        }
    }
    return ret;
}

template <typename KmerT>
void MatchFinder<KmerT>::setTiles(const flowcell::TileMetadataList &tiles)
{
    matchWriter_.reopen(iteration_, tiles);
}

/**
 * \brief Find matches.
 *
 * \param seeds Expects sorted ABCD permutation of seeds. Ns not included.
 */
template<typename KmerT>
const std::vector<MatchDistribution> &MatchFinder<KmerT>::findMatches(
    const typename std::vector<SeedT>::iterator seedsBegin,
    const std::vector<typename std::vector<SeedT>::iterator> &referenceSeedBounds,
    const bool findNeighbors,
    const bool finalPass)
{
    return  match(seedsBegin, referenceSeedBounds, findNeighbors, finalPass);
}

template<typename KmerT>
const std::vector<MatchDistribution> & MatchFinder<KmerT>::match(
    typename std::vector<SeedT>::const_iterator seedsBegin,
    const std::vector<typename std::vector<SeedT>::iterator> &referenceSeedBounds,
    const bool findNeighbors,
    const bool finalPass)
{

    // parallelize the matching across all the masks of all the references.
    typename KmerSourceMetadataList::const_iterator kmerSourceIterator = kmerSourceMetadataList_.begin();
    threads_.execute(boost::bind(&MatchFinder::matchMaskParallel, this,
                                 boost::ref(seedsBegin),
                                 boost::cref(referenceSeedBounds),
                                 findNeighbors,
                                 finalPass,
                                 boost::ref(kmerSourceIterator),
                                 _1),
                     threadsMax_);

    return threadMatchDistributions_;
}

template <typename KmerT>
std::pair<typename std::vector<typename MatchFinder<KmerT>::SeedT>::const_iterator, typename std::vector<typename MatchFinder<KmerT>::SeedT>::const_iterator>
MatchFinder<KmerT>::skipToTheNextMask(
    const typename std::vector<SeedT>::const_iterator currentBegin,
    const typename std::vector<SeedT>::const_iterator seedsEnd,
    const KmerT currentMask,
    const unsigned maskWidth,
    const bool storeNSeedNoMatches)
{
    const KmerT endSeed =
        (currentMask << (oligo::KmerTraits<KmerT>::KMER_BITS - maskWidth)) |
        ((~KmerT(0) << (oligo::KmerTraits<KmerT>::KMER_BITS - maskWidth)) ^ ~KmerT(0));
    // the first possibly-N-seed of the current mask (only the 111.. masks will actually have N-seeds, but it does not matter)
    const SeedT firstPossiblyNSeed(endSeed, SMALLEST_N_SEED_ID);
    typename std::vector<SeedT>::const_iterator currentEnd =
        std::lower_bound(currentBegin, seedsEnd, firstPossiblyNSeed, &orderByKmerSeedIndex<KmerT>);
    // roll to the first non-N seed so that next thread does not have to deal with our Ns
    typename std::vector<SeedT>::const_iterator nextBegin =
        std::find_if(currentEnd, seedsEnd, boost::bind(&SeedT::isNSeed, _1) != true);

    if (nextBegin != currentEnd)
    {
        ISAAC_THREAD_CERR << "Skipped " << nextBegin - currentEnd << " N-seeds for mask " << boost::numeric_cast<int>(currentMask) << std::endl;
        if (storeNSeedNoMatches)
        {
            // generate no-match entries for N-containing seeds or else the match selector statistics will never see those clusters
            BOOST_FOREACH(const SeedT &seed, std::make_pair(currentEnd, nextBegin))
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

template <typename KmerT>
void MatchFinder<KmerT>::matchMaskParallel(
    typename std::vector<SeedT>::const_iterator &seedsBegin,
    const std::vector<typename std::vector<SeedT>::iterator> &referenceSeedBounds,
    const bool findNeighbors,
    const bool finalPass,
    typename KmerSourceMetadataList::const_iterator &kmerSourceIterator,
    const unsigned threadNumber)
{
    boost::lock_guard<boost::mutex> lock(mutex_);
    while (kmerSourceIterator < kmerSourceMetadataList_.end() && referenceSeedBounds.back() != seedsBegin)
    {
        typename KmerSourceMetadataList::const_iterator ourKmerSource = kmerSourceIterator++;
        const KmerT currentMask = ourKmerSource->mask_;

        const typename std::vector<SeedT>::const_iterator ourBegin(seedsBegin);
        const typename std::vector<SeedT>::const_iterator seedsEnd(referenceSeedBounds.at(ourKmerSource->referenceIndex_));

        // on the final pass make sure the n-seeds of the open reads get their
        // nomatches stored. Else Match selector stats will report incorrect total cluster count
        std::pair<typename std::vector<SeedT>::const_iterator, typename std::vector<SeedT>::const_iterator> ourEndNextBegin =
            skipToTheNextMask(ourBegin, seedsEnd, currentMask, ourKmerSource->maskWidth_, finalPass);
        seedsBegin = ourEndNextBegin.second;

        const boost::filesystem::path &sortedReferencePath = ourKmerSource->maskFilePath_;
        std::istream threadReferenceFile(threadReferenceFileBuffers_.at(threadNumber).get(sortedReferencePath, io::FileBufWithReopen::SequentialOften));

        {
            common::unlock_guard<boost::mutex> unlock(mutex_);

            if (findNeighbors)
            {
                matchFinder::NeighborMaskMatcher<KmerT>(
                    ignoreRepeats_,
                    repeatThreshold_,  neighborhoodSizeThreshold_,
                    seedMetadataList_,
                    referenceContigKaryotypes_.at(ourKmerSource->referenceIndex_),
                    foundExactMatchesOnly_).matchNeighborsMask(
                        ourBegin, ourEndNextBegin.first, currentMask,
                        threadMatchDistributions_[threadNumber],
                        threadRepeatLists_[threadNumber],
                        threadNeighborsLists_[threadNumber],
                        matchWriter_,
                        threadReferenceFile);
            }
            else
            {
                matchFinder::ExactMaskMatcher<KmerT>(
                    1 == iteration_, finalPass,
                    repeatThreshold_,
                    ignoreNeighbors_, seedMetadataList_,
                    referenceContigKaryotypes_.at(ourKmerSource->referenceIndex_),
                    foundExactMatchesOnly_).matchMask(
                        ourBegin, ourEndNextBegin.first, currentMask,
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

template class MatchFinder<oligo::ShortKmerType>;
template class MatchFinder<oligo::KmerType>;
template class MatchFinder<oligo::LongKmerType>;

} // namespace alignment
} // namespace isaac
