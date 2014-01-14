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
 ** \file FragmentBuilder.cpp
 **
 ** \brief See FragmentBuilder.hh
 ** 
 ** \author Come Raczy
 **/

#include <numeric>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "alignment/Alignment.hh"
#include "alignment/FragmentBuilder.hh"
#include "alignment/Quality.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

FragmentBuilder::FragmentBuilder(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const unsigned repeatThreshold,
    const unsigned maxSeedsPerRead,
    const unsigned gappedMismatchesMax,
    const bool avoidSmithWaterman,
    const int gapMatchScore,
    const int gapMismatchScore,
    const int gapOpenScore,
    const int gapExtendScore,
    const int minGapExtendScore,
    const unsigned gapLimit)
    : repeatThreshold_(repeatThreshold)
    , semialignedGapLimit_(gapLimit)
    , gappedMismatchesMax_(gappedMismatchesMax)
    , seedMatchCounts_(maxSeedsPerRead * readsMax_)
    , repeatSeedsCount_(0)
    , fragments_(readsMax_) // max number of read ends we ever have to deal with
    , cigarBuffer_(Cigar::getMaxOperationsForReads(flowcellLayoutList) *
                   // one seed generates up to repeat threshold matches for each strand
                   repeatThreshold_ * maxSeedsPerRead * 2)
    , ungappedAligner_(gapMatchScore, gapMismatchScore, gapOpenScore, gapExtendScore, minGapExtendScore)
    , gappedAligner_(flowcellLayoutList, avoidSmithWaterman, gapMatchScore, gapMismatchScore, gapOpenScore, gapExtendScore, minGapExtendScore)
    , simpleIndelAligner_(gapMatchScore, gapMismatchScore, gapOpenScore, gapExtendScore, minGapExtendScore, semialignedGapLimit_)
{
    std::for_each(fragments_.begin(), fragments_.end(),
                  boost::bind(&std::vector<FragmentMetadata>::reserve, _1,
                              // one seed generates up to repeat threshold matches for each strand
                              repeatThreshold_ * maxSeedsPerRead * 2));
}

void FragmentBuilder::clear()
{
    BOOST_FOREACH(std::vector<FragmentMetadata> &fragmentList, fragments_)
    {
        fragmentList.clear();
    }
    cigarBuffer_.clear();
    std::fill(seedMatchCounts_.begin(), seedMatchCounts_.end(), 0);
    repeatSeedsCount_ = 0;
}

/**
 * \return true if at least one fragment was built.
 */
bool FragmentBuilder::build(
    const std::vector<reference::Contig> &contigList,
    const flowcell::ReadMetadataList &readMetadataList,
    const SeedMetadataList &seedMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    std::vector<Match>::const_iterator matchBegin,
    const std::vector<Match>::const_iterator matchEnd,
    const Cluster &cluster,
    const bool withGaps)
{
    clear();
    if (matchBegin < matchEnd)
    {
        ISAAC_ASSERT_MSG(!matchBegin->isNoMatch(), "Fake match lists must be dealt with outside");
        ISAAC_ASSERT_MSG(cluster.getNonEmptyReadsCount() == readMetadataList.size(), "cluster geometry must match");
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "FragmentBuilder::build: cluster " << matchBegin->seedId.getCluster() << " (" << cluster.getId() << ")");

        for(;matchEnd != matchBegin && !matchBegin->isNoMatch(); ++matchBegin)
        {
            if (repeatThreshold_ > seedMatchCounts_.at(matchBegin->getSeedId().getSeed()))
            {
                if (matchBegin->isTooManyMatch())
                {
                    seedMatchCounts_[matchBegin->getSeedId().getSeed()] = repeatThreshold_;
                    ++repeatSeedsCount_;
                    const unsigned matchReadIndex = seedMetadataList[matchBegin->getSeedId().getSeed()].getReadIndex();
                    ISAAC_ASSERT_MSG(fragments_[matchReadIndex].empty(), "Too-many matches are expected to sort to the top");
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "FragmentBuilder::build: blocking read " << matchReadIndex  << ", seed " << matchBegin->getSeedId().getSeed());

                }
                else
                {
                    if (repeatThreshold_ == ++seedMatchCounts_[matchBegin->getSeedId().getSeed()])
                    {
                        ++repeatSeedsCount_;
                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "FragmentBuilder::build: read " << seedMetadataList[matchBegin->getSeedId().getSeed()].getReadIndex()  <<
                                                    ", seed " << matchBegin->getSeedId().getSeed() << " exceeded repeat threshold and will be ignored");
                    }
                    else
                    {
                        addMatch(readMetadataList, seedMetadataList, *matchBegin, cluster);
                    }
                }
            }
        }

        if (repeatSeedsCount_)
        {
            // if some of the seeds produced more than repeatThreshold_ matches,
            // we need to remove those matches to avoid bias
            std::for_each(fragments_.begin(), fragments_.end(),
                          boost::bind(&FragmentBuilder::removeRepeatSeedAlignments, this, _1));
        }

        if (fragments_.end() != std::find_if(fragments_.begin(), fragments_.end(),
                                             !boost::bind(&std::vector<FragmentMetadata>::empty, _1)))
        {
            alignFragments(contigList, readMetadataList, seedMetadataList, sequencingAdapters, withGaps);

            return true;
        }
    }
    return false;
}

void FragmentBuilder::alignFragments(
    const std::vector<reference::Contig> &contigList,
    const flowcell::ReadMetadataList &readMetadataList,
    const SeedMetadataList &seedMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const bool withGaps)
{
    unsigned fragmentIndex = 0;
    BOOST_FOREACH(std::vector<FragmentMetadata> &fragmentList, fragments_)
    {
        if (!fragmentList.empty())
        {
            consolidateDuplicateFragments(fragmentList, false);
            ISAAC_DEV_TRACE_BLOCK(const Cluster & cluster = fragmentList.at(0).getCluster();)
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "Ungapped alignment for cluster " << cluster.getTile() <<
                ":" << cluster.getId() << " [" << fragmentIndex << "] " << sequencingAdapters.size() << " adapters");
                
            matchSelector::FragmentSequencingAdapterClipper adapterClipper(sequencingAdapters);
            BOOST_FOREACH(FragmentMetadata &fragmentMetadata, fragmentList)
            {
                fragmentMetadata.repeatSeedsCount = repeatSeedsCount_;
                ISAAC_THREAD_CERR_DEV_TRACE("    Original   : " << fragmentMetadata);
                const unsigned contigId = fragmentMetadata.contigId;
                ISAAC_ASSERT_MSG(contigList.size() > contigId, "Unexpected contig id");
                const reference::Contig &contig = contigList[contigId];

                adapterClipper.checkInitStrand(fragmentMetadata, contig);
                ungappedAligner_.alignUngapped(fragmentMetadata, cigarBuffer_, readMetadataList, adapterClipper, contig);
                ISAAC_THREAD_CERR_DEV_TRACE("    Aligned    : " << fragmentMetadata);
            }

            // make sure all positions that are in the list are unique.
            consolidateDuplicateFragments(fragmentList, true);

            if (semialignedGapLimit_)
            {
                simpleIndelAligner_.alignSimpleIndels(cigarBuffer_, contigList, readMetadataList, seedMetadataList, fragmentList);
                consolidateDuplicateFragments(fragmentList, true);
            }

            ISAAC_THREAD_CERR_DEV_TRACE("Gapped alignment for cluster " << cluster.getTile() <<
                ":" << cluster.getId() << " [" << fragmentIndex << "] " << sequencingAdapters.size() << " adapters");
            // If there are still bad alignments, try to do expensive smith-waterman on them.
            BOOST_FOREACH(FragmentMetadata &fragmentMetadata, fragmentList)
            {
                const unsigned contigId = fragmentMetadata.contigId;
                const reference::Contig &contig = contigList[contigId];

                adapterClipper.checkInitStrand(fragmentMetadata, contig);
                ISAAC_THREAD_CERR_DEV_TRACE("    Original    : " << fragmentMetadata);
                if (withGaps && BandedSmithWaterman::mismatchesCutoff < fragmentMetadata.mismatchCount)
                {
                    FragmentMetadata tmp = fragmentMetadata;
                    const unsigned matchCount = gappedAligner_.alignGapped(tmp, cigarBuffer_, readMetadataList, adapterClipper, contig);
                    ISAAC_THREAD_CERR_DEV_TRACE("    Gap-aligned: " << tmp);
                    if (matchCount && matchCount + BandedSmithWaterman::WIDEST_GAP_SIZE > fragmentMetadata.getObservedLength() &&
                        (tmp.mismatchCount <= gappedMismatchesMax_) &&
                        (fragmentMetadata.mismatchCount > tmp.mismatchCount) &&
                        ISAAC_LP_LESS(fragmentMetadata.logProbability, tmp.logProbability))
                    {
                        ISAAC_THREAD_CERR_DEV_TRACE("    Using gap-aligned: " << tmp);
                        fragmentMetadata = tmp;
                    }
                }
            }
            // gapped alignment and adapter trimming may adjust the alignment position
            consolidateDuplicateFragments(fragmentList, true);
        }
        ++fragmentIndex;
    }
}

void FragmentBuilder::addMatch(
    const flowcell::ReadMetadataList &readMetadataList,
    const SeedMetadataList &seedMetadataList, const Match &match,
    const Cluster &cluster)
{
    const unsigned seedIndex = match.seedId.getSeed();
    const SeedMetadata &seedMetadata = seedMetadataList.at(seedIndex);
    const unsigned readIndex = seedMetadata.getReadIndex();
    const reference::ReferencePosition &seedLocation = match.location;
    const unsigned contigId = seedLocation.getContigId();
    const bool reverse = match.seedId.isReverse();
    const long readPosition = getReadPosition(readMetadataList, seedMetadata, seedLocation.getPosition(), reverse);

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "    adding match at position " << contigId << ":" << readPosition);
    // Cigar buffer and offset are set separately at the point when we align
    fragments_[readIndex].push_back(FragmentMetadata(&cluster, NULL, readIndex));
    FragmentMetadata &fragment = fragments_[readIndex].back();
    fragment.firstSeedIndex = seedIndex;
    fragment.contigId = contigId;
    fragment.position = readPosition;
    fragment.reverse = reverse;
    if (seedMetadata.getLength() != STRONG_SEED_LENGTH && seedLocation.hasNeighbors())
    {
        fragment.nonUniqueSeedOffsets.first = std::min<unsigned>(fragment.nonUniqueSeedOffsets.first, seedMetadata.getOffset());
        fragment.nonUniqueSeedOffsets.second = std::max<unsigned>(fragment.nonUniqueSeedOffsets.second, seedMetadata.getOffset());
    }
    else
    {
        fragment.uniqueSeedCount = 1;
    }
}


inline bool exceeedsThreshold(const std::vector<unsigned> &vect, const unsigned index, const unsigned threshold)
{
    return vect[index] >= threshold;
}
/**
 * \brief Removes entries introduced by the seeds which later on have been found to have over repeatThreshold_
 *        matches
 *
 */
void FragmentBuilder::removeRepeatSeedAlignments(FragmentMetadataList &fragmentList)
{
    fragmentList.erase(std::remove_if(fragmentList.begin(), fragmentList.end(),
                                      boost::bind(&exceeedsThreshold, boost::ref(seedMatchCounts_), boost::bind(&FragmentMetadata::getFirstSeedIndex, _1), repeatThreshold_)),
                                      fragmentList.end());
}

/**
 * \brief Removes entries designating the same alignment location and strand
 *
 * This is important not only to avoid repetitive processing but also to avoid good alignments
 * scored poorly because they are present multiple times in the list.
 *
 * \param removeUnaligned if true, entries with !isAlignment() are removed
 *
 * \postcondition the fragmentList is ordered
 * \postcondition The fragmentList contains unique alignments only
 */
void FragmentBuilder::consolidateDuplicateFragments(FragmentMetadataList &fragmentList, const bool removeUnaligned)
{
    ISAAC_THREAD_CERR_DEV_TRACE("consolidateDuplicateFragments initial size: " << fragmentList.size());
    // although initially matches arrive ordered by location, gapped alignment might have moved the
    // start position of some.
    std::sort(fragmentList.begin(), fragmentList.end());
    std::vector<FragmentMetadata>::iterator lastFragment = fragmentList.begin();
    while(fragmentList.end() != lastFragment && removeUnaligned && !lastFragment->isAligned())
    {
        ++lastFragment;
    }

    lastFragment = fragmentList.erase(fragmentList.begin(), lastFragment);

    if (2 > fragmentList.size())
    {
        return;
    }

    ISAAC_ASSERT_MSG(fragmentList.end() != lastFragment, "Unexpected end in list of size " << fragmentList.size());
    ISAAC_ASSERT_MSG(fragmentList.end() != lastFragment + 1, "Unexpected end - 1 in list of size " << fragmentList.size());

    for (std::vector<FragmentMetadata>::iterator currentFragment = lastFragment + 1; fragmentList.end() != currentFragment; ++currentFragment)
    {
        if (removeUnaligned && !currentFragment->isAligned())
        {
            ISAAC_THREAD_CERR_DEV_TRACE("    ignoring " << *currentFragment);
        }
        else if (*lastFragment == *currentFragment)
        {
            ISAAC_THREAD_CERR_DEV_TRACE("    updating " << *lastFragment << " with " << *currentFragment);
            lastFragment->consolidate(*currentFragment);
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE("    adding   " << *currentFragment << " prev:" << *lastFragment);
            ++lastFragment;
            if (lastFragment != currentFragment)
            {
                *lastFragment = *currentFragment;
            }
        }
    }
    fragmentList.resize(1 + lastFragment - fragmentList.begin());
    ISAAC_THREAD_CERR_DEV_TRACE("consolidateDuplicateFragments consolidated size: " << fragmentList.size());
}

inline long FragmentBuilder::getReadPosition(
    const flowcell::ReadMetadataList &readMetadataList,
    const SeedMetadata &seedMetadata, const long seedPosition, const bool reverse) const
{
    const int seedOffset = seedMetadata.getOffset();
    if (reverse)
    {
        const unsigned readIndex = seedMetadata.getReadIndex();
        const unsigned readLength = readMetadataList.at(readIndex).getLength();
        const unsigned seedLength = seedMetadata.getLength();
        // 'seedPosition + seedLength + seedOffset' is the first position past the end of the read
        return seedPosition + seedLength + seedOffset - readLength;
    }
    else
    {
        return seedPosition - seedOffset;
    }
}

} // namespace alignment
} // namespace isaac
