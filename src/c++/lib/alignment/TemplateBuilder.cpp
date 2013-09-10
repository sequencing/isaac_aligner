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
 ** \file TemplateBuilder.cpp
 **
 ** \brief See TemplateBuilder.hh
 **
 ** \author Come Raczy
 **/

#include <numeric>
#include <limits>

#include <boost/function_output_iterator.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "alignment/Alignment.hh"
#include "alignment/FragmentBuilder.hh"
#include "alignment/Quality.hh"
#include "alignment/TemplateBuilder.hh"
#include "oligo/KmerGenerator.hpp"

namespace isaac
{
namespace alignment
{

const unsigned TemplateBuilder::DODGY_BUT_CLEAN_ALIGNMENT_SCORE;

/**
 * \brief checks if the alignment simply does not make sense regardless
 *        of how unique it is.
 *
 * A bad alignment is either:
 * - less than 32 consecutive matches and either
 * - an alignment with too many mismatches (more than 1/8 the of the unclipped length) or
 * - an alignment with a low log probability (arbitrarily set to the equivalent of a read with
 *   all bases being Q40 and mismatches on 1/4 of them)
 */
inline bool isVeryBadAlignment(const FragmentMetadata &fragment)
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "fragment.mismatchCount > fragment.getMappedLength() / 8 " << fragment.mismatchCount << ">" << fragment.getMappedLength() / 8);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "fragment.logProbability < LOG_MISMATCH_Q40 / 4 * fragment.getMappedLength() " << fragment.logProbability << "<" << LOG_MISMATCH_Q40 / 4 * fragment.getMappedLength());
    return fragment.matchesInARow < 32 &&
        (fragment.mismatchCount > fragment.getMappedLength() / 8 || fragment.logProbability < LOG_MISMATCH_Q40 / 4 * fragment.getMappedLength());
}

TemplateBuilder::TemplateBuilder(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const unsigned repeatThreshold,
    const unsigned maxSeedsPerRead,
    const bool scatterRepeats,
    const unsigned gappedMismatchesMax,
    const bool avoidSmithWaterman,
    const int gapMatchScore,
    const int gapMismatchScore,
    const int gapOpenScore,
    const int gapExtendScore,
    const int minGapExtendScore,
    const unsigned semialignedGapLimit,
    const DodgyAlignmentScore dodgyAlignmentScore)
    : scatterRepeats_(scatterRepeats)
    , dodgyAlignmentScore_(dodgyAlignmentScore)
    , fragmentBuilder_(flowcellLayoutList, repeatThreshold, maxSeedsPerRead, gappedMismatchesMax,
                       avoidSmithWaterman, gapMatchScore, gapMismatchScore, gapOpenScore, gapExtendScore, minGapExtendScore, semialignedGapLimit)
    , bamTemplate_(fragmentBuilder_.getCigarBuffer())
    , shadowAligner_(flowcellLayoutList,
                     gappedMismatchesMax, avoidSmithWaterman, gapMatchScore, gapMismatchScore, gapOpenScore, gapExtendScore, minGapExtendScore)
    , cigarBuffer_(10000)
    , shadowList_(TRACKED_REPEATS_MAX_ONE_READ)
{
    // In the weirdest case each seed will generate --repeat-threshold alignments which then will generate TRACKED_REPEATS_MAX_ONE_READ shadows
    // So, max seeds per read times max number of orphans times the max number of shadows that can be rescued for each orphan times the number of reads
    allPairProbabilities_.reserve(maxSeedsPerRead * repeatThreshold * TRACKED_REPEATS_MAX_ONE_READ * readsMax_);
    for (unsigned r = 0; r < readsMax_; ++r)
    {
        // In the weirdest case each seed will generate --repeat-threshold alignments which then will generate TRACKED_REPEATS_MAX_ONE_READ shadows
        // So, max seeds per read times the max number of orphans times the max number of shadows that can be rescued for each orphan
        // plus the max number of shadows that have been discovered during seed matching...
        allShadowProbabilities_[r].reserve(maxSeedsPerRead *
                                           repeatThreshold * TRACKED_REPEATS_MAX_ONE_READ +
                                           TRACKED_REPEATS_MAX_ONE_READ);
    }
}
bool TemplateBuilder::buildTemplate(
    const std::vector<reference::Contig> &contigList,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const Cluster &cluster,
    const TemplateLengthStatistics &templateLengthStatistics,
    const unsigned mapqThreshold)
{
    const std::vector<std::vector<FragmentMetadata> > &fragments = fragmentBuilder_.getFragments();
    bool ret = buildTemplate(
        contigList, restOfGenomeCorrection, readMetadataList, sequencingAdapters, fragments, cluster, templateLengthStatistics);

    if (ret && bamTemplate_.hasAlignmentScore())
    {
        if (!bamTemplate_.isProperPair())
        {
            // for improper pairs, mark fragment individually unaligned if they are below threshold
            ret = bamTemplate_.filterLowQualityFragments(mapqThreshold);
        }
        else if (mapqThreshold > bamTemplate_.getAlignmentScore())
        {
            // mark the whole thing unaligned if the pair is below threshold
            bamTemplate_.filterLowQualityFragments(-1U);
            ret = false;
        }
    }
    return ret;
}
bool TemplateBuilder::buildTemplate(
    const std::vector<reference::Contig> &contigList,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    const Cluster &cluster,
    const TemplateLengthStatistics &templateLengthStatistics)
{
    cigarBuffer_.clear();
    // Initialize the bamTemplate_ with unaligned fragments for the cluster
    bamTemplate_.initialize(readMetadataList, cluster);
    if (2 == readMetadataList.size() && 2 == fragments.size())
    {
        if (!fragments[0].empty() && !fragments[1].empty())
        {
            return pickBestPair(
                contigList, restOfGenomeCorrection, readMetadataList, sequencingAdapters, fragments, templateLengthStatistics);
        }
        else if (!fragments[0].empty() || !fragments[1].empty())
        {
            return rescueShadow(
                contigList, restOfGenomeCorrection, readMetadataList, sequencingAdapters, fragments, templateLengthStatistics);
        }
        else
        {
            // TODO: implement a recovery mechanism for unaligned clusters
            return false;
        }
    }
    else if(1 == readMetadataList.size() || 1 == fragments.size())
    {
        ISAAC_ASSERT_MSG(fragments[1].empty(), "With single-ended data expecting the fragment to be placed at index 0");
        if (!fragments[0].empty())
        {
            return pickBestFragment(restOfGenomeCorrection, templateLengthStatistics, fragments[0]);
        }
        else
        {
            // TODO: implement a recovery mechanism for unaligned clusters
            return false;
        }
    }
    else
    {
        ISAAC_ASSERT_MSG(false, (boost::format("TemplateBuilder supports at most 2 reads: %d reads found") % fragments.size()).str().c_str());
        return false;
    }
}


std::vector<FragmentMetadata>::const_iterator TemplateBuilder::getBestFragment(const std::vector<FragmentMetadata> &fragmentList) const
{

/*
    assert(!fragmentList.empty());
    using boost::lambda::bind;
    using boost::lambda::_1;
    using boost::lambda::_2;
    return std::min_element(
        fragmentList.begin(),
        fragmentList.end(),
        bind<double>(&FragmentMetadata::smithWatermanScore, _1) < bind<double>(&FragmentMetadata::smithWatermanScore, _2));

*/


    common::FiniteCapacityVector<std::vector<FragmentMetadata>::const_iterator, TRACKED_REPEATS_MAX_ONE_READ> bestFragments;

    unsigned bestFragmentScore = -1U;
    double bestFragmentLogProbability = -std::numeric_limits<double>::max();
    for(FragmentIterator fragmentIterator = fragmentList.begin();
        fragmentList.end() != fragmentIterator; ++fragmentIterator)
    {
        if (bestFragmentScore > fragmentIterator->smithWatermanScore ||
            (bestFragmentScore == fragmentIterator->smithWatermanScore &&
                ISAAC_LP_LESS(bestFragmentLogProbability, fragmentIterator->logProbability)))
        {
            bestFragmentScore = fragmentIterator->smithWatermanScore;
            bestFragmentLogProbability = fragmentIterator->logProbability;
            bestFragments.clear();
            bestFragments.push_back(fragmentIterator);
        }
        else if (bestFragmentScore == fragmentIterator->smithWatermanScore &&
            ISAAC_LP_EQUALS(bestFragmentLogProbability, fragmentIterator->logProbability))
        {
            bestFragments.push_back(fragmentIterator);
        }

    }

    const unsigned clusterId = fragmentList[0].getCluster().getId();
    const unsigned repeatIndex = scatterRepeats_ ? (clusterId % bestFragments.size()) : 0;

    ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::getBestFragment returning repeat " <<
                                repeatIndex << " out of " << bestFragments.size() << " cluster id: " << clusterId << " " <<
                                *bestFragments[repeatIndex]);

    return bestFragments[repeatIndex];

}

/**
 * \brief Set mapping score given the probability of the best choice and probabilities of alternative choices
 *
 * \return true if the end is considered well anchored.
 */
bool TemplateBuilder::updateMappingScore(
    FragmentMetadata &fragment,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const TemplateLengthStatistics &templateLengthStatistics,
    const FragmentIterator listFragment,
    const std::vector<FragmentMetadata> &fragmentList,
    const bool forceWellAnchored) const
{
    ISAAC_ASSERT_MSG(fragmentList.end() == std::adjacent_find(fragmentList.begin(), fragmentList.end()), "Expecting unique alignments in fragment list");
    if (forceWellAnchored || fragment.isWellAnchored())
    {
        // Either seed(s) without neighbors or mutliple neighbor-having seeds agree
        // NOTE: the assumption is that if two (or more) separate seeds designate the same alignment position,
        // their neighbors are very unlikely to do so
        double neighborProbability = restOfGenomeCorrection.getReadRogCorrection(listFragment->getReadIndex());
        unsigned long bestScore = 0;
        unsigned bestHitCount = 0;
        unsigned otherHitsCount = 0;
        for (std::vector<FragmentMetadata>::const_iterator i = fragmentList.begin(); fragmentList.end() != i; ++i)
        {
            if (!bestHitCount || i->smithWatermanScore < bestScore)
            {
                otherHitsCount += bestHitCount;
                bestHitCount = 1;
                bestScore = i->smithWatermanScore;
            }
            else if (i->smithWatermanScore == bestScore)
            {
                ++bestHitCount;
            }
            else
            {
                ++otherHitsCount;
            }

            if (listFragment != i)
            {
                neighborProbability += exp(i->logProbability);
            }
        }
        fragment.alignmentScore = unsigned(floor(-10.0 *  log10(neighborProbability / (neighborProbability + exp(listFragment->logProbability)))));
        ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::updateMappingScore: " << bestHitCount << ":" << otherHitsCount <<
                                    " sm=" << fragment.getAlignmentScore());
        return true;
    }
    else
    {
        // All seeds have neighbors and no consensus among them.
        fragment.alignmentScore = 0;
        ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::updateMappingScore: all seeds have neighbors and neither two agree");
        return false;
    }
}

void TemplateBuilder::locateBestPair(
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    const TemplateLengthStatistics &templateLengthStatistics,
    TemplateBuilder::BestPairInfo &ret) const
{
    ret.init(fragments[0].begin(), fragments[1].begin());
    typedef std::vector<FragmentMetadata>::const_iterator FragmentIterator;
    FragmentIterator contigBegin[2] = {fragments[0].begin(), fragments[1].begin()};
    FragmentIterator contigEnd[2];
    // process the pairs of fragments contigId by contigId
    while ((fragments[0].end() != contigBegin[0]) && (fragments[1].end() != contigBegin[1]))
    {
        for (std::size_t i = 0; 2 > i; ++i)
        {
            contigEnd[i] = contigBegin[i] + 1;
            while((fragments[i].end() != contigEnd[i]) && (contigEnd[i]->contigId == contigBegin[i]->contigId))
            {
                ++contigEnd[i];
            }
            if (fragments[i].end() != contigEnd[i])
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(677606/*fragments[0].front().getCluster().getId()*/, "locateBestPair: read " << i << ": " << *contigBegin[i] << "-"  << *contigEnd[i]);
            }
            else
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), "locateBestPair: read " << i << ": " << *contigBegin[i] << "-END");
            }
        }
        if (contigBegin[0]->contigId == contigBegin[1]->contigId)
        {
            // pick best pair for the contig
            ISAAC_THREAD_CERR_DEV_TRACE((boost::format("locateBestPair: contig %d") % contigBegin[1]->contigId).str());
            FragmentIterator currentFragment[2] = {contigBegin[0], contigBegin[1]};
            while (contigEnd[0] != currentFragment[0])
            {
                currentFragment[1] = contigBegin[1];
                while(contigEnd[1] != currentFragment[1])
                {
                    assert(currentFragment[0]->contigId == currentFragment[1]->contigId);
                    if (templateLengthStatistics.matchModel(*currentFragment[0], *currentFragment[1]))
                    {
                        const double currentLogProbability = currentFragment[0]->logProbability + currentFragment[1]->logProbability;
                        const double currentProbability = exp(currentLogProbability);
                        const unsigned long templateScore = currentFragment[0]->smithWatermanScore + currentFragment[1]->smithWatermanScore;
                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), "locateBestPair pair: " <<
                            *currentFragment[0] << "-" << *currentFragment[1] << " (" <<
                            currentFragment[0]->logProbability + currentFragment[1]->logProbability << "lp. " << templateScore << "ubs)");
                        ret.totalTemplateProbability += currentProbability;
                        if (0 == ret.resolvedTemplateCount ||
                            ret.bestTemplateScore > templateScore ||
                            (templateScore == ret.bestTemplateScore && ISAAC_LP_LESS(ret.bestTemplateLogProbability, currentLogProbability)))
                        {
                            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), "reset");
                            ret.bestPairFragments[0].clear();
                            ret.bestPairFragments[1].clear();
                            ret.bestPairFragments[0].push_back(currentFragment[0]);
                            ret.bestPairFragments[1].push_back(currentFragment[1]);

                            ret.bestTemplateScore = templateScore;
                            ret.bestTemplateLogProbability = currentLogProbability;

                        }
                        else if (templateScore == ret.bestTemplateScore && ISAAC_LP_EQUALS(currentLogProbability, ret.bestTemplateLogProbability))
                        {
                            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), "append");
                            // this will blow up if we go over TRACKED_REPEATS_MAX_EVER, but at the moment
                            // it is unclear what to do best in such situation. Deal with it when it happens...
                            ret.bestPairFragments[0].push_back(currentFragment[0]);
                            ret.bestPairFragments[1].push_back(currentFragment[1]);
                        }

                        ++ret.resolvedTemplateCount;


                    }
                    ++currentFragment[1];
                }

                ++currentFragment[0];
            }
            // advance both lists to the next contig
            for (std::size_t i = 0; 2 > i; ++i)
            {
                contigBegin[i] = contigEnd[i];
            }
        }
        else
        {
            // advance the list with the smallest contigId to the next contig
            const std::size_t i = contigBegin[0]->contigId < contigBegin[1]->contigId ? 0 : 1;
            contigBegin[i] = contigEnd[i];
        }
    }
    if (ret.resolvedTemplateCount)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(
            fragments[0].front().getCluster().getId(),
            "locateBestPair " << ret << " blp=" << ret.bestPairFragments[0][0]->logProbability + ret.bestPairFragments[1][0]->logProbability);
        ret.bestPairEditDistance = ret.bestPairFragments[0][0]->getEditDistance() + ret.bestPairFragments[1][0]->getEditDistance();
    }
    else
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), "locateBestPair: nothing good...");
    }
}

/**
 * Computes alignment scores. If multiple best pair combinations are possible, and repeat scattering is
 * enabled, brings the selected repeat on top of the bestPairInfo.bestPairFragments
 * \return  true if the template is well anchored and no realignment seems to be needed.
 */
bool TemplateBuilder::buildPairedEndTemplate(
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const TemplateLengthStatistics &templateLengthStatistics,
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    BestPairInfo &bestPairInfo)
{
    // initialize the bamTemplate_ with the selected fragments
    FragmentMetadata &read1 = bamTemplate_.getFragmentMetadata(0);
    FragmentMetadata &read2 = bamTemplate_.getFragmentMetadata(1);

    if (scatterRepeats_)
    {
        // pick a pseudo-random template out of the list of equivalently probable ones
        const unsigned repeatIndex = read1.getCluster().getId() % bestPairInfo.bestPairFragments[0].size();
        ISAAC_ASSERT_MSG(bestPairInfo.bestPairFragments[0].size() == bestPairInfo.bestPairFragments[1].size(),
                         "Each pair of fragments with the same index must represent a single template");
        // put it on top
        std::swap(bestPairInfo.bestPairFragments[0][0], bestPairInfo.bestPairFragments[0][repeatIndex]);
        std::swap(bestPairInfo.bestPairFragments[1][0], bestPairInfo.bestPairFragments[1][repeatIndex]);

        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(read1.getCluster().getId(), "TemplateBuilder::buildPairedEndTemplate: Picked repeat " << repeatIndex <<
                                    " out of " << bestPairInfo.bestPairFragments[0].size() << " " << bamTemplate_);
    }

    read1 = *bestPairInfo.bestPairFragments[0][0];
    read2 = *bestPairInfo.bestPairFragments[1][0];

    const bool r1WellAnchored =
        updateMappingScore(read1, restOfGenomeCorrection, templateLengthStatistics, bestPairInfo.bestPairFragments[0][0], fragments[0], read2.isWellAnchored());
    const bool r2WellAnchored =
        updateMappingScore(read2, restOfGenomeCorrection,  templateLengthStatistics, bestPairInfo.bestPairFragments[1][0], fragments[1], read1.isWellAnchored());

    bamTemplate_.setProperPair(alignment::TemplateLengthStatistics::Nominal == templateLengthStatistics.checkModel(read1, read2));
    if (r1WellAnchored || r2WellAnchored)
    {
        const double otherPairsProbability =
            (bestPairInfo.totalTemplateProbability - exp(bestPairInfo.bestTemplateLogProbability)) +
            restOfGenomeCorrection.getRogCorrection();
        bamTemplate_.setAlignmentScore(
            unsigned(floor(-10.0 * log10(otherPairsProbability /
                                         (bestPairInfo.totalTemplateProbability +
                                             restOfGenomeCorrection.getRogCorrection())))));

        ISAAC_ASSERT_MSG(bestPairInfo.bestPairFragments[0].size() == 1 || bamTemplate_.getAlignmentScore() < 4,
                         (boost::format("alignment score too high for a repeat of %d: %s") % bestPairInfo.bestPairFragments[0].size() % bamTemplate_).str().c_str());

        if (r1WellAnchored && r2WellAnchored &&
            !read1.repeatSeedsCount && !read2.repeatSeedsCount)
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(read1.getCluster().getId(), "Pair-end  template well anchored: " << bamTemplate_);
            // Pair sticks quite well, no need to try rescuing.
            return true;
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(read1.getCluster().getId(), "Pair-end  template only one end is well anchored: " << bamTemplate_);
            // One or both of the ends dodgy. See what rescuing gives
            return false;
        }
    }
    else
    {
        bamTemplate_.setAlignmentScore(-1U);
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(read1.getCluster().getId(), "Pair-end  template looks quite random: " << bamTemplate_);
        return false;
    }

}

bool TemplateBuilder::flagDodgyTemplate(FragmentMetadata &orphan, FragmentMetadata &shadow, BamTemplate &bamTemplate) const
{
    if (DODGY_ALIGNMENT_SCORE_UNALIGNED == dodgyAlignmentScore_)
    {
        // both must sort into unaligned bin. setUnaligned will not do it.
        orphan.setNoMatch();
        shadow.setNoMatch();
        bamTemplate.setAlignmentScore(-1U);
        return false;
    }
    else if (DODGY_ALIGNMENT_SCORE_UNKNOWN == dodgyAlignmentScore_)
    {
        orphan.alignmentScore = -1U;
        shadow.alignmentScore = -1U;
        bamTemplate.setAlignmentScore(-1U);
    }
    else
    {
    	// numeric score assigned. Keep it at -1U now. Bam generator will figure out.
    	orphan.alignmentScore = -1U;
    	shadow.alignmentScore = -1U;
    	bamTemplate.setAlignmentScore(-1U);
    }

    ISAAC_THREAD_CERR_DEV_TRACE("flagDodgyTemplate: " << bamTemplate_);
    return true;
}

bool TemplateBuilder::rescueShadow(
    const std::vector<reference::Contig> &contigList,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    const TemplateLengthStatistics &templateLengthStatistics)
{
    const unsigned orphanIndex = fragments[0].empty() ? 1 : 0;
    const unsigned shadowIndex = (orphanIndex + 1) % readsMax_;
    typedef std::vector<FragmentMetadata>::const_iterator FragmentIterator;
    const FragmentIterator bestOrphanIterator = getBestFragment(fragments[orphanIndex]);
    BestPairInfo &bestPair = bestRescuedPair_;
    bestPair.clear();
    bestPair.bestPairFragments[orphanIndex].push_back(bestOrphanIterator);

    allShadowProbabilities_[orphanIndex].clear();
    for(FragmentIterator orphanIterator = fragments[orphanIndex].begin();
        fragments[orphanIndex].end() != orphanIterator; ++orphanIterator)
    {
        const FragmentMetadata &orphan = *orphanIterator;
        //const double orphanProbability = exp(orphan.logProbability);
        shadowList_.clear();
        ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::rescueShadow Orphan: " << orphan);

        if (ISAAC_LP_LESS(orphan.logProbability + orphanLogProbabilitySlack_, bestOrphanIterator->logProbability))
        {
            ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::rescueShadow orphan too bad to try rescuing shadows");
        }
        else if (shadowAligner_.rescueShadow(contigList, orphan, shadowList_, readMetadataList, sequencingAdapters,
                                        templateLengthStatistics, 0))
        {
            // the best shadow for this orphan is the first in the list
            const FragmentMetadata &bestRescued = shadowList_.front();

            const double currentTemplateLogProbability = orphan.logProbability + bestRescued.logProbability;

            const unsigned long templateScore = orphan.smithWatermanScore + bestRescued.smithWatermanScore;

            if(isVeryBadAlignment(bestRescued))
            {
                ISAAC_THREAD_CERR_DEV_TRACE("    Rescued shadow too bad: " << bestRescued);
            }
            else
            {
                if (0 == bestPair.resolvedTemplateCount || templateScore < bestPair.bestTemplateScore ||
                    (templateScore == bestPair.bestTemplateScore && ISAAC_LP_LESS(bestPair.bestTemplateLogProbability, currentTemplateLogProbability)))
                {
                    bestPair.bestTemplateLogProbability = currentTemplateLogProbability;
                    bestPair.bestTemplateScore = templateScore;
                    bestPair.bestPairFragments[orphanIndex].clear();
                    bestPair.bestPairFragments[orphanIndex].push_back(orphanIterator);

                    bestOrphanShadows_[orphanIndex].clear();
                    bestOrphanShadows_[orphanIndex].push_back(cloneWithCigar(bestRescued));

                }
                else if (templateScore == bestPair.bestTemplateScore && ISAAC_LP_EQUALS(currentTemplateLogProbability, bestPair.bestTemplateLogProbability))
                {
                    bestPair.bestPairFragments[orphanIndex].push_back(orphanIterator);
                    bestOrphanShadows_[orphanIndex].push_back(cloneWithCigar(bestRescued));
                }

                ++bestPair.resolvedTemplateCount;
            }

        }
        else if (!shadowList_.empty())
        {
            // The shadow hits a repetitive region next to one of the orphans
            ISAAC_THREAD_CERR_DEV_TRACE((boost::format("Shadow rescue hits a repeat. Orphan: %s") % orphan).str());
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE("No shadows rescued");
        }

        BOOST_FOREACH(const FragmentMetadata &shadow, shadowList_)
        {
            allShadowProbabilities_[orphanIndex].push_back(ShadowProbability(shadow));
            bestPair.totalTemplateProbability += exp(orphan.logProbability + shadow.logProbability);
        }

    }

    const double totalShadowProbability = (0 < bestPair.resolvedTemplateCount) ?
        sumUniqueShadowProbabilities(allShadowProbabilities_[orphanIndex].begin(), allShadowProbabilities_[orphanIndex].end()) : 0.0;
    
    bool ret = true;
    FragmentMetadata &orphan = bamTemplate_.getFragmentMetadata(orphanIndex);
    if (0 < bestPair.resolvedTemplateCount)
    {
        const unsigned clusterId = fragments[orphanIndex][0].getCluster().getId();
        const unsigned repeatIndex = scatterRepeats_ ? clusterId % bestPair.bestPairFragments[orphanIndex].size() : 0;

        orphan = *bestPair.bestPairFragments[orphanIndex][repeatIndex];
        FragmentMetadata &bestShadow = bestOrphanShadows_[orphanIndex][repeatIndex];

        const bool assumeWellAnchored = updateMappingScore(
            orphan, restOfGenomeCorrection, templateLengthStatistics, bestPair.bestPairFragments[orphanIndex][repeatIndex],
            fragments[orphanIndex],
            // if pair perfectly matches give it a chance to have some alignment score
            0 == orphan.getEditDistance() + bestShadow.getEditDistance());
        if (assumeWellAnchored)
        {
            const double shadowRog = restOfGenomeCorrection.getReadRogCorrection(bestShadow.getReadIndex());

            // if the orphan aligned in a reasonably unique way, use probabilities to compute alignment scores
            const double otherShadowsProbability = (totalShadowProbability - exp(bestShadow.logProbability)) + shadowRog;
            bestShadow.alignmentScore = floor(-10.0 * log10(otherShadowsProbability / (totalShadowProbability + shadowRog)));

            const double otherPairsProbability = (bestPair.totalTemplateProbability - exp(bestPair.bestTemplateLogProbability)) +
                restOfGenomeCorrection.getRogCorrection();
            bamTemplate_.setAlignmentScore(unsigned(floor(-10.0 * log10(otherPairsProbability / (bestPair.totalTemplateProbability +
                restOfGenomeCorrection.getRogCorrection())))));

            if (!orphan.alignmentScore || !orphan.isWellAnchored())
            {
                // CASAVA variant caller is sensitive to coverage dips. For pairs that are not
                // anchored well but don't produce any variants, ensure the alignment score does not get too high
                // but does not cause CASAVA to ignore their existence (default alignment score cutoff in starling is 10)
                bamTemplate_.setAlignmentScore(std::min(DODGY_BUT_CLEAN_ALIGNMENT_SCORE, bamTemplate_.getAlignmentScore()));
                bestShadow.alignmentScore = std::min(DODGY_BUT_CLEAN_ALIGNMENT_SCORE, bestShadow.alignmentScore);
                orphan.alignmentScore = std::min(DODGY_BUT_CLEAN_ALIGNMENT_SCORE, orphan.alignmentScore);
            }

            ISAAC_ASSERT_MSG(bestPair.bestPairFragments[orphanIndex].size() == 1 || bamTemplate_.getAlignmentScore() < 4,
                             (boost::format("alignment score too high for a repeat of %d: %s") % bestPair.bestPairFragments[orphanIndex].size() % bamTemplate_).str().c_str());
        }
        else
        {
            ret = flagDodgyTemplate(orphan, bestShadow, bamTemplate_);
        }
        bamTemplate_.getFragmentMetadata(shadowIndex) = bestShadow;
        bamTemplate_.setProperPair(alignment::TemplateLengthStatistics::Nominal == templateLengthStatistics.checkModel(orphan, bestShadow));
        ISAAC_THREAD_CERR_DEV_TRACE("rescueShadow: rescued  template: " << bamTemplate_);
        if (scatterRepeats_)
        {
            ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::rescueShadow: Picked repeat " << repeatIndex <<
                                        " out of " << bestPair.bestPairFragments[orphanIndex].size() << " " << bamTemplate_);
        }

    }
    else
    {
        orphan = *bestOrphanIterator;
        FragmentMetadata &shadow = bamTemplate_.getFragmentMetadata(shadowIndex);

        if (isVeryBadAlignment(orphan))
        {
            ISAAC_THREAD_CERR_DEV_TRACE("    Singleton too bad: " << orphan);
            // both must sort into unaligned bin. setUnaligned will not do it.
            orphan.setNoMatch();
            shadow.setNoMatch();
            ret = false;
        }
        else
        {
            // mark shadow as 'shadow' (singleton's position, etc)
            shadow.contigId = orphan.contigId;
            shadow.position = orphan.position;
            shadow.readIndex = shadowIndex;
            shadow.alignmentScore = 0;
            shadow.cigarLength = 0;
            if (!updateMappingScore(orphan, restOfGenomeCorrection, templateLengthStatistics, bestOrphanIterator, fragments[orphanIndex], 0 == orphan.getEditDistance()))
            {
                ret = flagDodgyTemplate(orphan, shadow, bamTemplate_);
            }
            else
            {
                if (!orphan.isWellAnchored())
                {
                    orphan.setAlignmentScore(std::min(DODGY_BUT_CLEAN_ALIGNMENT_SCORE, orphan.getAlignmentScore()));
                }
                bamTemplate_.setAlignmentScore(0);
                ISAAC_THREAD_CERR_DEV_TRACE("rescueShadow: keeping unrescued: " << bamTemplate_);
            }
        }
    }

    return ret;
}

FragmentMetadata TemplateBuilder::cloneWithCigar(const FragmentMetadata &right)
{
    FragmentMetadata ret = right;
    ret.cigarBuffer = &cigarBuffer_;
    ret.cigarOffset = cigarBuffer_.size();
    const std::vector<unsigned>::const_iterator cigarBegin = right.cigarBuffer->begin() + right.cigarOffset;
    const std::vector<unsigned>::const_iterator cigarEnd = cigarBegin + right.cigarLength;
    cigarBuffer_.insert(cigarBuffer_.end(), cigarBegin, cigarEnd);
    return ret;
}

inline void addToDouble(double &to, double what)
{
    to += what;
}

double TemplateBuilder::sumUniqueShadowProbabilities(std::vector<ShadowProbability>::iterator begin, std::vector<ShadowProbability>::iterator end)
{
    double ret = 0.0;
    std::sort(begin, end);

    std::unique_copy(begin, end, boost::make_function_output_iterator(
        bind(&addToDouble, boost::ref(ret), boost::bind(&exp, boost::bind(&ShadowProbability::logProbability, _1)))));

    return ret;
}

double TemplateBuilder::sumUniquePairProbabilities(std::vector<PairProbability>::iterator begin, std::vector<PairProbability>::iterator end)
{
    double ret = 0.0;
    std::sort(begin, end);

    std::unique_copy(begin, end, boost::make_function_output_iterator(
        bind(&addToDouble, boost::ref(ret), boost::bind(&exp, boost::bind(&PairProbability::logProbability, _1)))));

    return ret;
}

bool TemplateBuilder::buildDisjoinedTemplate(
    const std::vector<reference::Contig> &contigList,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    const TemplateLengthStatistics &templateLengthStatistics,
    const BestPairInfo &knownBestPair)
{
    const FragmentIterator bestDisjoinedFragments[readsMax_] = {getBestFragment(fragments[0]), getBestFragment(fragments[1])};
    const long bestTemplateLength = knownBestPair.getBestTemplateLength();

    unsigned bestOrphanIndex = 0;
    // bestOrphans contains an iterator to the read used for template anchoring on each of the read numbers.
    // bestOrphanIndex tells which read number anchors the best one. Note that if there is a number of equivalently
    // good templates (in case of repeat) bestOrphanIndex contains all of them even if they are originally anchored
    // by different read number reads
    BestPairInfo &bestOrphans = bestRescuedPair_;
    bestOrphans.init(bestDisjoinedFragments[0], bestDisjoinedFragments[1]);

    allPairProbabilities_.clear();
    for (unsigned orphanIndex = 0; readsMax_ > orphanIndex; ++orphanIndex)
    {
        allShadowProbabilities_[orphanIndex].clear();
        bestOrphanShadows_[orphanIndex].clear();

        for(FragmentIterator orphanIterator = fragments[orphanIndex].begin();
            fragments[orphanIndex].end() != orphanIterator;
            ++orphanIterator)
        {
            const FragmentMetadata &orphan = *orphanIterator;

            const bool skipThisOrphan = (knownBestPair.resolvedTemplateCount ? //when ed cutoff is set, no point in trying orphans with greater edit distance
                    orphan.getEditDistance() > (knownBestPair.bestPairEditDistance + SKIP_ORPHAN_EDIT_DISTANCE) :
                    ISAAC_LP_LESS(orphan.logProbability + orphanLogProbabilitySlack_, bestDisjoinedFragments[orphanIndex]->logProbability));

            shadowList_.clear();
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "TemplateBuilder::buildDisjoinedTemplate Orphan: " << orphan);
            if (!skipThisOrphan && shadowAligner_.rescueShadow(contigList, orphan, shadowList_, readMetadataList,
                                                               sequencingAdapters,
                                                               templateLengthStatistics, bestTemplateLength))
            {
                const FragmentMetadata &bestRescued = shadowList_.front();
                const double currentTemplateLogProbability = orphan.logProbability + bestRescued.logProbability;

                const unsigned rescuedEditDistance = orphan.getEditDistance() + bestRescued.getEditDistance();

                if (isVeryBadAlignment(bestRescued))
                {
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescued shadow too bad: " << bestRescued);
                }
                else if (!knownBestPair.resolvedTemplateCount || (knownBestPair.bestPairEditDistance + SKIP_ORPHAN_EDIT_DISTANCE) >= rescuedEditDistance)
                {
                    const unsigned long templateScore = orphan.smithWatermanScore + bestRescued.smithWatermanScore;

                    if (0 == bestOrphans.resolvedTemplateCount ||
                        templateScore < bestOrphans.bestTemplateScore ||
                        (templateScore == bestOrphans.bestTemplateScore && ISAAC_LP_LESS(bestOrphans.bestTemplateLogProbability, currentTemplateLogProbability)))
                    {
                        bestOrphans.bestTemplateLogProbability = currentTemplateLogProbability;
                        bestOrphans.bestTemplateScore = templateScore;

                        bestOrphans.bestPairFragments[orphanIndex].clear();

                        bestOrphans.bestPairFragments[orphanIndex].push_back(orphanIterator);
                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "new: " << orphanIndex <<
                                                    "oi " << templateScore << "ts " << bestOrphans.bestTemplateScore << "bts " <<
                                                    currentTemplateLogProbability << "ctlp " << bestOrphans.bestTemplateLogProbability << "btlp");


                        // bestOrphanShadows_ contains the shadows that correspond to each  of bestOrphans[bestOrphanIndex]
                        bestOrphanShadows_[orphanIndex].clear();
                        bestOrphanShadows_[orphanIndex].push_back(cloneWithCigar(bestRescued));

                        bestOrphanIndex = orphanIndex;
                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "buildDisjoinedTemplate: candidate template: " << *orphanIterator <<
                                                    "-" << bestRescued <<
                                                    " for edit distance " << rescuedEditDistance <<
                                                    " <= " << knownBestPair.bestPairEditDistance);
                    }
                    else if (templateScore == bestOrphans.bestTemplateScore && ISAAC_LP_EQUALS(currentTemplateLogProbability, bestOrphans.bestTemplateLogProbability))
                    {
                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "append: " << orphanIndex <<
                                                    "oi " << templateScore << "ts " << bestOrphans.bestTemplateScore << "bts " <<
                                                    currentTemplateLogProbability << "ctlp " << bestOrphans.bestTemplateLogProbability << "btlp " <<
                                                    bestOrphans.bestPairFragments[orphanIndex].size() << "bpfsize");
                        bestOrphans.bestPairFragments[orphanIndex].push_back(orphanIterator);
                        bestOrphanShadows_[orphanIndex].push_back(cloneWithCigar(bestRescued));
                    }
                    ++bestOrphans.resolvedTemplateCount;
                }
                else
                {
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "buildDisjoinedTemplate: ignoring template: " << orphan <<
                        "-" << bestRescued << " for edit distance " << rescuedEditDistance);
                }
            }
            else if (!shadowList_.empty())
            {
                // The shadow hits a repetitive region next to one of the orphans
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "Mate rescue hits a repeat. Orphan: " << orphan);
            }

            BOOST_FOREACH(const FragmentMetadata &shadow, shadowList_)
            {
                allPairProbabilities_.push_back(PairProbability(0 == orphanIndex ? orphan : shadow, 0 == orphanIndex ? shadow : orphan));
                allShadowProbabilities_[orphanIndex].push_back(ShadowProbability(shadow));

/*
                ISAAC_THREAD_CERR_DEV_TRACE("apped probabilities " << (orphan.logProbability + shadow.logProbability) <<
                                            " " << (0 == orphanIndex ? orphan.getFStrandReferencePosition() : shadow.getFStrandReferencePosition()) <<
                                            " " << (0 == orphanIndex ? shadow.getFStrandReferencePosition() : orphan.getFStrandReferencePosition()));
*/

            }
        }
    }

    const unsigned bestShadowIndex = (bestOrphanIndex + 1) % readsMax_;
    double totalShadowProbability = 0.0;
    double totalOrphanProbability = 0.0;

    if (0 < bestOrphans.resolvedTemplateCount)
    {
        // mix in all the seed-discovered shadows before computing the totalShadowProbability
        BOOST_FOREACH(const FragmentMetadata &shadow, fragments[bestShadowIndex])
        {
            allShadowProbabilities_[bestOrphanIndex].push_back(ShadowProbability(shadow));
        }
        totalShadowProbability = sumUniqueShadowProbabilities(allShadowProbabilities_[bestOrphanIndex].begin(), allShadowProbabilities_[bestOrphanIndex].end());

        BOOST_FOREACH(const FragmentMetadata &orphan, fragments[bestOrphanIndex])
        {
            allShadowProbabilities_[bestShadowIndex].push_back(ShadowProbability(orphan));
        }
        totalOrphanProbability = sumUniqueShadowProbabilities(allShadowProbabilities_[bestShadowIndex].begin(), allShadowProbabilities_[bestShadowIndex].end());

        bestOrphans.totalTemplateProbability += sumUniquePairProbabilities(allPairProbabilities_.begin(), allPairProbabilities_.end());
    }

    return scoreDisjoinedTemplate(
        fragments,
        restOfGenomeCorrection,
        templateLengthStatistics,
        bestOrphans,
        knownBestPair,
        bestOrphanIndex,
        totalShadowProbability,
        totalOrphanProbability,
        bestDisjoinedFragments);
}

bool TemplateBuilder::scoreDisjoinedTemplate(
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const TemplateLengthStatistics &templateLengthStatistics,
    const BestPairInfo &bestOrphans,
    const BestPairInfo &knownBestPair,
    const unsigned bestOrphanIndex,
    const double totalShadowProbability,
    const double totalOrphanProbability,
    const FragmentIterator bestDisjoinedFragments[readsMax_])
{
    bool ret = true;
    if (0 < bestOrphans.resolvedTemplateCount)
    {
        const unsigned clusterId = fragments[0][0].getCluster().getId();
        const unsigned repeatIndex = scatterRepeats_ ? clusterId % bestOrphans.bestPairFragments[bestOrphanIndex].size() : 0;

        const FragmentMetadata &bestOrphan = *bestOrphans.bestPairFragments[bestOrphanIndex][repeatIndex];
        FragmentMetadata &bestShadow = bestOrphanShadows_[bestOrphanIndex][repeatIndex];

        // rediscovered is most likely going to be false for repeat pairs, which is ok as we don't want repeats to get high scores anyway
        const bool rediscovered = !repeatIndex && knownBestPair.resolvedTemplateCount &&
            *knownBestPair.bestPairFragments[bestOrphan.getReadIndex()][0] == bestOrphan  &&
            *knownBestPair.bestPairFragments[bestShadow.getReadIndex()][0] == bestShadow;

        FragmentMetadata &orphan = bamTemplate_.getFragmentMetadata(bestOrphan.getReadIndex());
        orphan = bestOrphan;
        const bool shadowWellAnchored = rediscovered && knownBestPair.bestPairFragments[bestShadow.getReadIndex()][0]->isWellAnchored();
        const bool assumeWellAnchored =
            updateMappingScore(
                orphan, restOfGenomeCorrection, templateLengthStatistics,
                bestOrphans.bestPairFragments[bestOrphan.getReadIndex()][repeatIndex],
                fragments[bestOrphan.getReadIndex()],
                // allow to have alignment score if the whole pair
                // does not produce false variants or the shadow was well anchored
                0 == orphan.getEditDistance() + bestShadow.getEditDistance() || shadowWellAnchored);

        bamTemplate_.setProperPair(alignment::TemplateLengthStatistics::Nominal == templateLengthStatistics.checkModel(orphan, bestShadow));
        if (assumeWellAnchored)
        {
            const double shadowRog = restOfGenomeCorrection.getReadRogCorrection(bestShadow.getReadIndex());
            // orphan aligns uniquely enough to use probabilities for shadow and template alignment scores
            // notice the subtraction ensured to happen before adding shadowRog as shadowRog tends to be very small
            const double otherShadowsProbability = (totalShadowProbability - exp(bestShadow.logProbability)) + shadowRog;
            bestShadow.alignmentScore = floor(-10.0 * log10(otherShadowsProbability / (totalShadowProbability + shadowRog)));

            const double orphanRog = restOfGenomeCorrection.getReadRogCorrection(bestOrphan.getReadIndex());
            const double otherOrphansProbability = (totalOrphanProbability - exp(bestOrphan.logProbability)) + orphanRog;
            orphan.alignmentScore = floor(-10.0 * log10(otherOrphansProbability / (totalOrphanProbability + orphanRog)));

            const double otherPairsProbability = (bestOrphans.totalTemplateProbability - exp(bestOrphans.bestTemplateLogProbability)) +
                restOfGenomeCorrection.getRogCorrection();
            bamTemplate_.setAlignmentScore(unsigned(floor(-10.0 * log10(otherPairsProbability / (bestOrphans.totalTemplateProbability +
                restOfGenomeCorrection.getRogCorrection())))));

            if ((!orphan.alignmentScore || !orphan.isWellAnchored()) &&
                (!bestShadow.alignmentScore || !shadowWellAnchored))
            {
                ISAAC_THREAD_CERR_DEV_TRACE("buildDisjoinedTemplate: rescuing unanchored template: " << bamTemplate_ << bestOrphans.bestTemplateLogProbability << ":" << bestOrphans.totalTemplateProbability << "blp:tp");
                // CASAVA variant caller is sensitive to coverage dips. For pairs that are not
                // anchored well but don't produce any variants, ensure the alignment score does not get too high
                // but does not cause CASAVA to ignore their existence (default alignment score cutoff in starling is 10)
                bamTemplate_.setAlignmentScore(std::min(DODGY_BUT_CLEAN_ALIGNMENT_SCORE, bamTemplate_.getAlignmentScore()));
                bestShadow.alignmentScore = std::min(DODGY_BUT_CLEAN_ALIGNMENT_SCORE, bestShadow.alignmentScore);
                orphan.alignmentScore = std::min(DODGY_BUT_CLEAN_ALIGNMENT_SCORE, orphan.alignmentScore);
                bamTemplate_.getFragmentMetadata(bestShadow.getReadIndex()) = bestShadow;
                ISAAC_THREAD_CERR_DEV_TRACE("buildDisjoinedTemplate: rescued  unanchored template: " << bamTemplate_ << bestOrphans.bestTemplateLogProbability << ":" << bestOrphans.totalTemplateProbability << "blp:tp");
            }
            else
            {
                bamTemplate_.getFragmentMetadata(bestShadow.getReadIndex()) = bestShadow;
                ISAAC_THREAD_CERR_DEV_TRACE("buildDisjoinedTemplate: rescued  anchored template: " <<
                                            bamTemplate_ << bestOrphans.bestTemplateLogProbability << ":" << bestOrphans.totalTemplateProbability << "blp:tp " <<
                                            bestOrphans.bestPairFragments[bestOrphanIndex].size() << "bpfsize");
            }

            ISAAC_ASSERT_MSG(bestOrphans.bestPairFragments[bestOrphanIndex].size() == 1 || bamTemplate_.getAlignmentScore() < 4,
                             (boost::format("alignment score too high for a repeat of %d: %s") %
                                 bestOrphans.bestPairFragments[bestOrphanIndex].size() % bamTemplate_).str().c_str() <<
                                 *bestOrphans.bestPairFragments[bestOrphanIndex][repeatIndex] <<
                                 std::string(bestOrphans.bestPairFragments[bestOrphanIndex][repeatIndex]->getRead().getForwardSequence().begin(),
                                             bestOrphans.bestPairFragments[bestOrphanIndex][repeatIndex]->getRead().getForwardSequence().end()) << "," <<
                                 bestOrphanShadows_[bestOrphanIndex][repeatIndex] <<
                                 std::string(bestOrphanShadows_[bestOrphanIndex][repeatIndex].getRead().getForwardSequence().begin(),
                                             bestOrphanShadows_[bestOrphanIndex][repeatIndex].getRead().getForwardSequence().end()));
        }
        else
        {
            ret = flagDodgyTemplate(orphan, bestShadow, bamTemplate_);
            bamTemplate_.getFragmentMetadata(bestShadow.getReadIndex()) = bestShadow;
            ISAAC_THREAD_CERR_DEV_TRACE("buildDisjoinedTemplate: rescued  template: " << bamTemplate_);
        }

        if (scatterRepeats_)
        {
            ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::buildDisjoinedTemplate: Picked repeat " << repeatIndex <<
                                        " out of " << bestOrphans.bestPairFragments[bestOrphanIndex].size() << " " << bamTemplate_);
        }
    }
    else if (knownBestPair.resolvedTemplateCount)
    {
        // don't build disjoined template if there is a reasonable pair with knownBestPair.bestPairEditDistance already
        // However, since we have not been able to rediscover this pair, it means it is either outside of
        // consensus template boundaries or keeps hitting highly repetitive  locations with the shadow  all the time.
        ISAAC_THREAD_CERR_DEV_TRACE("Disjoined template: nothing better but the rediscovery did not happen");
        ret = flagDodgyTemplate(bamTemplate_.getFragmentMetadata(0), bamTemplate_.getFragmentMetadata(1), bamTemplate_);
    }
    else
    {
        FragmentMetadata &read1 = bamTemplate_.getFragmentMetadata(0);
        FragmentMetadata &read2 = bamTemplate_.getFragmentMetadata(1);
        read1 = *bestDisjoinedFragments[0];
        read2 = *bestDisjoinedFragments[1];
        bamTemplate_.setAlignmentScore(0);
        bamTemplate_.setProperPair(false);

        const bool assumeR1WellAnchored =
            updateMappingScore(read1, restOfGenomeCorrection, templateLengthStatistics, bestDisjoinedFragments[0], fragments[0], 0 == read1.getEditDistance());
        const bool assumeR2WellAnchored =
            updateMappingScore(read2, restOfGenomeCorrection, templateLengthStatistics, bestDisjoinedFragments[1], fragments[1], 0 == read2.getEditDistance());

        if (!assumeR1WellAnchored && !assumeR2WellAnchored)
        {
            // For anomalous pairs keep only if both ends is anchored well.
            ret = flagDodgyTemplate(read1, read2, bamTemplate_);
        }
        else
        {
            if (!read1.isWellAnchored())
            {
                read1.setAlignmentScore(std::min(DODGY_BUT_CLEAN_ALIGNMENT_SCORE, read1.getAlignmentScore()));
            }
            if (!read2.isWellAnchored())
            {
                read2.setAlignmentScore(std::min(DODGY_BUT_CLEAN_ALIGNMENT_SCORE, read2.getAlignmentScore()));
            }
            ISAAC_THREAD_CERR_DEV_TRACE("Disjoined template: keeping disjoined: " << bamTemplate_);
        }
    }
    return ret;
}

bool TemplateBuilder::flagDodgyTemplate(FragmentMetadata &orphan, BamTemplate &bamTemplate) const
{
    if (DODGY_ALIGNMENT_SCORE_UNALIGNED == dodgyAlignmentScore_)
    {
        // both must sort into unaligned bin. setUnaligned will not do it.
        orphan.setNoMatch();
        bamTemplate.setAlignmentScore(-1U);
        return false;
    }
    else if (DODGY_ALIGNMENT_SCORE_UNKNOWN == dodgyAlignmentScore_)
    {
        orphan.alignmentScore = -1U;
        bamTemplate.setAlignmentScore(-1U);
    }
    else
    {
        // numeric score assigned. Keep it at -1U now. Bam generator will figure out.
        orphan.alignmentScore = -1U;
        bamTemplate.setAlignmentScore(-1U);
    }

    ISAAC_THREAD_CERR_DEV_TRACE("flagDodgyTemplate: " << bamTemplate_);
    return true;
}

bool TemplateBuilder::pickBestFragment(
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const TemplateLengthStatistics &templateLengthStatistics,
    const std::vector<FragmentMetadata> &fragmentList)
{
    if (!fragmentList.empty())
    {
        typedef std::vector<FragmentMetadata>::const_iterator FragmentIterator;
        const FragmentIterator bestFragment = getBestFragment(fragmentList);
        bamTemplate_.getFragmentMetadata(0) = *bestFragment;
        if (!updateMappingScore(bamTemplate_.getFragmentMetadata(0), restOfGenomeCorrection, templateLengthStatistics, bestFragment, fragmentList, false))
        {
            return flagDodgyTemplate(bamTemplate_.getFragmentMetadata(0), bamTemplate_);
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate_.getFragmentMetadata(0).getCluster().getId(), "Orphan template: " << bamTemplate_.getFragmentMetadata(0));
            return true;
        }
    }

    ISAAC_THREAD_CERR_DEV_TRACE("Orphan template not aligned: " << bamTemplate_.getFragmentMetadata(0));
    return false;
}

bool TemplateBuilder::pickBestPair(
    const std::vector<reference::Contig> &contigList,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    const TemplateLengthStatistics &templateLengthStatistics)
{
    ISAAC_ASSERT_MSG(readsMax_ == fragments.size(), "TemplateBuilder::pickBestPair must be called for paired templates only");
    ISAAC_ASSERT_MSG(!fragments[0].empty() && !fragments[1].empty(), "TemplateBuilder::pickBestPair must be called for paired templates only");

    locateBestPair(fragments, templateLengthStatistics, bestCombinationPairInfo_);

    if (!bestCombinationPairInfo_.resolvedTemplateCount ||
        !buildPairedEndTemplate(restOfGenomeCorrection, templateLengthStatistics, fragments, bestCombinationPairInfo_) ||
        bestCombinationPairInfo_.bestPairEditDistance)
    {
        // nothing resolved from match combinations or the resolved pair does not have enough
        // unique matches to anchor it in a trustworthy way. Give rescuing a chance to find something
        // better or lower-scored
        return buildDisjoinedTemplate(contigList, restOfGenomeCorrection, readMetadataList, sequencingAdapters, fragments,
                               templateLengthStatistics, bestCombinationPairInfo_);
    }

    return true;

}

} // namespace alignment
} // namespace isaac
