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
 ** \file MatchSelector.cpp
 **
 ** Component to select the best matches among all possible candidates.
 **
 ** \author Come Raczy
 **/

#include <numeric>
#include <fstream>
#include <cerrno>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/thread.hpp>

#include "alignment/Alignment.hh"
#include "alignment/MatchSelector.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/FastIo.hh"
#include "reference/Contig.hh"
#include "reference/ContigLoader.hh"

#include "alignment/matchSelector/MatchSelectorStatsXml.hh"

namespace isaac
{
namespace alignment
{

bool orderByTotalReadLengthDesc(const flowcell::FlowcellLayoutList &flowcellLayoutList,
                                const flowcell::TileMetadata &left, const flowcell::TileMetadata &right)
{
    const unsigned leftTotalReadLength = flowcell::getTotalReadLength(flowcellLayoutList.at(left.getFlowcellIndex()).getReadMetadataList());
    const unsigned rightTotalReadLength = flowcell::getTotalReadLength(flowcellLayoutList.at(right.getFlowcellIndex()).getReadMetadataList());
    // Also keep the natural order of tiles when the read lengths are the same so that it is easier to track the
    // progress by monitoring the log output
    // This is important for fastq as there we don't have freedom to change the order of tile processing within one lane
    ISAAC_ASSERT_MSG(left.getFlowcellId() != right.getFlowcellId() || leftTotalReadLength == rightTotalReadLength,
                     "Tiles of the same flowcell must have the same read lengths: " << left << " vs " << right);
    return  leftTotalReadLength > rightTotalReadLength ||
        (leftTotalReadLength == rightTotalReadLength && left.getIndex() < right.getIndex());
}

const flowcell::TileMetadataList sortByTotalReadLengthDesc(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    flowcell::TileMetadataList tileMetadataList)
{
    std::sort(tileMetadataList.begin(), tileMetadataList.end(),
              boost::bind(&orderByTotalReadLengthDesc, boost::ref(flowcellLayoutList), _1, _2));
    return tileMetadataList;
}

std::vector<matchSelector::SequencingAdapterList> generateSequencingAdapters(const flowcell::BarcodeMetadataList &barcodeMetadataList)
{
    std::vector<matchSelector::SequencingAdapterList> ret(barcodeMetadataList.size());

    std::vector<matchSelector::SequencingAdapterList>::iterator barcodeAdaptersIterator = ret.begin();
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        BOOST_FOREACH(const flowcell::SequencingAdapterMetadata &adapter, barcode.getAdapters())
        {
            barcodeAdaptersIterator->push_back(matchSelector::SequencingAdapter(adapter));
        }
        ++barcodeAdaptersIterator;
    }

    return ret;
}

struct MatchDistributionContigFilter
{
    const MatchDistribution &matchDistribution_;
    MatchDistributionContigFilter(const MatchDistribution &matchDistribution) : matchDistribution_(matchDistribution){}
    bool isMapped(unsigned, unsigned contigIndex) const {return !matchDistribution_.isEmptyContig(contigIndex);}
};

MatchSelector::MatchSelector(
        matchSelector::FragmentStorage &fragmentStorage,
        const MatchDistribution &matchDistribution,
        const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const unsigned int maxThreadCount,
        const TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned repeatThreshold,
        const int mateDriftRange,
        const TemplateLengthStatistics &userTemplateLengthStatistics,
        const unsigned mapqThreshold,
        const bool perTileTls,
        const bool pfOnly,
        const unsigned baseQualityCutoff,
        const bool keepUnaligned,
        const bool clipSemialigned,
        const bool clipOverlapping,
        const bool scatterRepeats,
        const unsigned gappedMismatchesMax,
        const bool avoidSmithWaterman,
        const int gapMatchScore,
        const int gapMismatchScore,
        const int gapOpenScore,
        const int gapExtendScore,
        const int minGapExtendScore,
        const unsigned semialignedGapLimit,
        const TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore
    )
    : computeThreads_(maxThreadCount),
      tileMetadataList_(tileMetadataList),
      barcodeMetadataList_(barcodeMetadataList),
      flowcellLayoutList_(flowcellLayoutList),
      repeatThreshold_(repeatThreshold),
      userTemplateLengthStatistics_(userTemplateLengthStatistics),
      mapqThreshold_(mapqThreshold),
      perTileTls_(perTileTls),
      pfOnly_(pfOnly),
      baseQualityCutoff_(baseQualityCutoff),
      keepUnaligned_(keepUnaligned),
      clipSemialigned_(clipSemialigned),
      clipOverlapping_(clipOverlapping),
      barcodeSequencingAdapters_(generateSequencingAdapters(barcodeMetadataList_)),
      allStats_(tileMetadataList_.size(), matchSelector::MatchSelectorStats(barcodeMetadataList_)),
      threadStats_(computeThreads_.size(), matchSelector::MatchSelectorStats(barcodeMetadataList_)),
      matchDistribution_(matchDistribution),
      contigList_(reference::loadContigs(sortedReferenceMetadataList, MatchDistributionContigFilter(matchDistribution_), computeThreads_)),
      fragmentStorage_(fragmentStorage),
      threadCluster_(computeThreads_.size(),
                     Cluster(flowcell::getMaxReadLength(flowcellLayoutList_) +
                             flowcell::getMaxBarcodeLength(flowcellLayoutList_))),
      threadTemplateBuilders_(computeThreads_.size()),
      threadSemialignedEndsClippers_(clipSemialigned_ ? computeThreads_.size() : 0),
      threadOverlappingEndsClippers_(computeThreads_.size()),
      templateLengthDistribution_(mateDriftRange)
{
    while(threadTemplateBuilders_.size() < computeThreads_.size())
    {
        threadTemplateBuilders_.push_back(new TemplateBuilder(flowcellLayoutList_,
                                                              repeatThreshold_,
                                                              flowcell::getMaxSeedsPerRead(flowcellLayoutList_),
                                                              scatterRepeats,
                                                              gappedMismatchesMax,
                                                              avoidSmithWaterman,
                                                              gapMatchScore,
                                                              gapMismatchScore,
                                                              gapOpenScore,
                                                              gapExtendScore,
                                                              minGapExtendScore,
                                                              semialignedGapLimit,
                                                              dodgyAlignmentScore));
    }

    templateLengthDistribution_.reserve(flowcell::getMaxTileClusters(tileMetadataList_));

    ISAAC_THREAD_CERR << "Constructed the match selector" << std::endl;
}

void MatchSelector::dumpStats(const boost::filesystem::path &statsXmlPath)
{
    std::for_each(allStats_.begin(), allStats_.end(), boost::bind(&matchSelector::MatchSelectorStats::finalize, _1));

    // xml tree serialization used to take quite a bit of ram (now it does not). make sure it's available anyway
    {
        std::vector<matchSelector::MatchSelectorStats>().swap(threadStats_);
    }

    std::ofstream os(statsXmlPath.string().c_str());
    if (!os) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: Unable to open file for writing: " + statsXmlPath.string()));
    }

    matchSelector::MatchSelectorStatsXml statsXml(flowcellLayoutList_, barcodeMetadataList_, tileMetadataList_, allStats_);
    statsXml.serialize(os);
}

TemplateLengthStatistics MatchSelector::determineTemplateLength(
    const flowcell::TileMetadata &tileMetadata,
    const std::vector<reference::Contig> &barcodeContigList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const std::vector<Match>::const_iterator barcodeMatchListBegin,
    const std::vector<Match>::const_iterator barcodeMatchListEnd,
    const BclClusters &bclData,
    const unsigned threadNumber)
{
    const flowcell::Layout &flowcell = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex());
    const flowcell::ReadMetadataList &tileReads = flowcell.getReadMetadataList();
    templateLengthDistribution_.reset(barcodeContigList, tileReads);

    ISAAC_ASSERT_MSG(2 >= tileReads.size(), "only single-ended and paired reads are supported");

    if (2 != tileReads.size())
    {
        ISAAC_THREAD_CERR << "Using unstable template-length statistics for single-ended data: " << templateLengthDistribution_.getStatistics() << std::endl;
        return templateLengthDistribution_.getStatistics();
    }

    if (userTemplateLengthStatistics_.isStable())
    {
        ISAAC_THREAD_CERR << "Using user-defined template-length statistics: " << userTemplateLengthStatistics_ << std::endl;
        return userTemplateLengthStatistics_;
    }
    else
    {
        const SeedMetadataList &tileSeeds = flowcell.getSeedMetadataList();
        const unsigned barcodeLength = flowcell.getBarcodeLength();
        TemplateBuilder &ourThreadTemplateBuilder = threadTemplateBuilders_.at(threadNumber);
        Cluster& ourThreadCluster = threadCluster_.at(threadNumber);

        for (std::vector<Match>::const_iterator matchBegin(barcodeMatchListBegin), matchEnd(findNextCluster(barcodeMatchListBegin, barcodeMatchListEnd));
            barcodeMatchListEnd != matchBegin && !templateLengthDistribution_.getStatistics().isStable();
            matchBegin = matchEnd, matchEnd = findNextCluster(matchBegin, barcodeMatchListEnd))
        {
            // identify all the matches for the current cluster
            const unsigned int clusterId = matchBegin->getCluster();
            ISAAC_ASSERT_MSG(clusterId < tileMetadata.getClusterCount(), "Cluster ids are expected to be 0-based within the tile.");

            // use only good pf clusters for template length calculation (and no fake matchlists)
            if (bclData.pf(clusterId) && !matchBegin->location.isNoMatch())
            {
                // initialize the cluster with the bcl data
                ourThreadCluster.init(tileReads, bclData.cluster(clusterId),
                                      matchBegin->getTile(), matchBegin->getCluster(),
                                      bclData.xy(clusterId), true, barcodeLength);
                // build the fragments for that cluster

/*
                // prevent clipping during template length statistics calculation
                static const matchSelector::SequencingAdapterList noSequencingAdapters;
                ourThreadTemplateBuilder.buildFragments(barcodeContigList, tileReads, tileSeeds, noSequencingAdapters,
                                                        matchBegin, matchEnd, ourThreadCluster, false);
*/
                ourThreadTemplateBuilder.buildFragments(barcodeContigList, tileReads, tileSeeds, sequencingAdapters,
                                                        matchBegin, matchEnd, ourThreadCluster, false);
                // use the fragments to build the template length statistics
                templateLengthDistribution_.addTemplate(ourThreadTemplateBuilder.getFragments());
            }
        }
        if (!templateLengthDistribution_.isStable())
        {
            templateLengthDistribution_.finalize();
        }
    }
    return templateLengthDistribution_.getStatistics();
}

void MatchSelector::processMatchList(
    const std::vector<reference::Contig> &barcodeContigList,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const std::pair<std::vector<Match>::const_iterator, std::vector<Match>::const_iterator> ourMatchListBeginEnd,
    const flowcell::TileMetadata & tileMetadata,
    const BclClusters &bclData,
    const TemplateLengthStatistics & templateLengthStatistics,
    const unsigned threadNumber)
{

    Cluster &ourThreadCluster = threadCluster_[threadNumber];
    TemplateBuilder &ourThreadTemplateBuilder = threadTemplateBuilders_.at(threadNumber);
    matchSelector::MatchSelectorStats &ourThreadStats = threadStats_.at(threadNumber);
    BamTemplate &ourThreadBamTemplate = ourThreadTemplateBuilder.getBamTemplate();

    const flowcell::Layout &flowcell = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex());
    const SeedMetadataList &tileSeeds = flowcell.getSeedMetadataList();
    const flowcell::ReadMetadataList &tileReads = flowcell.getReadMetadataList();
    const std::size_t barcodeLength = flowcell.getBarcodeLength();

    unsigned uniqueClustersToSkip = computeThreads_.size() - threadNumber - 1;
    unsigned clusterId = ourMatchListBeginEnd.first->getCluster();
    for (std::vector<Match>::const_iterator matchBegin = ourMatchListBeginEnd.first; ourMatchListBeginEnd.second != matchBegin; ++matchBegin)
    {
        // skip to the first match of our clusterId
        if (matchBegin->getCluster() != clusterId)
        {
            --uniqueClustersToSkip;
            clusterId = matchBegin->getCluster();
        }
        if (!uniqueClustersToSkip)
        {
            uniqueClustersToSkip = computeThreads_.size();

            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(matchBegin->seedId.getCluster(), "MatchSelector::processMatchList: cluster " << matchBegin->seedId.getCluster());

            ISAAC_ASSERT_MSG(clusterId < tileMetadata.getClusterCount(), "Cluster ids are expected to be 0-based within the tile.");

            // initialize the cluster with the bcl data
            ourThreadCluster.init(tileReads, bclData.cluster(clusterId), matchBegin->getTile(), clusterId,
                                  bclData.xy(clusterId), bclData.pf(clusterId), barcodeLength);
            trimLowQualityEnds(ourThreadCluster, baseQualityCutoff_);

            if ((pfOnly_ && !bclData.pf(clusterId)) || matchBegin->location.isNoMatch())
            {
                // if pfOnly_ is set, this non-pf cluster will not be reported as a regularly-processed one.
                // if match list begins with noMatchReferencePosition, then this cluster does not have any matches at all. This is
                // because noMatchReferencePosition has the highest possible contig number and sort will put it to the end of match list
                // In etiher case report it as skipped to ensure statistics consistency
                ourThreadBamTemplate.initialize(tileReads, ourThreadCluster);
                ourThreadStats.recordTemplate(tileReads, templateLengthStatistics, ourThreadBamTemplate, matchBegin->getBarcode(),
                                              matchBegin->location.isNoMatch() ?
                                                  matchBegin->getSeedId().isNSeedId() ?
                                                      matchSelector::Qc :
                                                      matchSelector::NmNm :
                                                  matchSelector::Filtered);
                if (keepUnaligned_ && (!pfOnly_ || bclData.pf(clusterId)))
                {
                    fragmentStorage_.add(ourThreadBamTemplate, matchBegin->getBarcode());
                }
            }
            else //if the cluster has matches and is not filtered-out by pf filtering
            {
                // find the first match that does not belong to our clusterId
                const std::vector<Match>::const_iterator matchEnd = findNextCluster(matchBegin, ourMatchListBeginEnd.second);

                // build the fragments for that cluster
                if (ourThreadTemplateBuilder.buildFragments(barcodeContigList, tileReads, tileSeeds, sequencingAdapters,
                                                            matchBegin, matchEnd, ourThreadCluster, true))
                {
                    ISAAC_ASSERT_MSG(2 >= ourThreadBamTemplate.getFragmentCount(), "only paired and singed ended data supported");

                    // build the template for the fragments
                    if (ourThreadTemplateBuilder.buildTemplate(
                        barcodeContigList, restOfGenomeCorrection, tileReads, sequencingAdapters,
                        ourThreadCluster, templateLengthStatistics,
                        mapqThreshold_) || keepUnaligned_)
                    {
                        if (clipSemialigned_)
                        {
                            threadSemialignedEndsClippers_[threadNumber].reset();
                            threadSemialignedEndsClippers_[threadNumber].clip(barcodeContigList, ourThreadBamTemplate);
                        }
                        if (clipOverlapping_)
                        {
                            threadOverlappingEndsClippers_[threadNumber].reset();
                            threadOverlappingEndsClippers_[threadNumber].clip(barcodeContigList, ourThreadBamTemplate);
                        }
                        fragmentStorage_.add(ourThreadBamTemplate, matchBegin->getBarcode());
                    }

                    ourThreadStats.recordTemplate(tileReads, templateLengthStatistics, ourThreadBamTemplate,
                                                  matchBegin->getBarcode(), matchSelector::Normal);
                }
                else
                {
                    ourThreadBamTemplate.initialize(tileReads, ourThreadCluster);
                    ourThreadStats.recordTemplate(tileReads, templateLengthStatistics, ourThreadBamTemplate,
                                                  matchBegin->getBarcode(), matchSelector::Rm);
                    if (keepUnaligned_)
                    {
                        fragmentStorage_.add(ourThreadBamTemplate, matchBegin->getBarcode());
                    }
                }
                // save some cluster counting
                matchBegin = matchEnd - 1;
            }
        }
    }
}

void MatchSelector::parallelSelect(
    const MatchTally &matchTally,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    const flowcell::TileMetadata &tileMetadata,
    std::vector<Match> &matchList,
    const BclClusters &bclData)
{
    std::for_each(threadStats_.begin(), threadStats_.end(), boost::bind(&matchSelector::MatchSelectorStats::reset, _1));

    ISAAC_THREAD_CERR << "Resizing fragment storage for " <<  tileMetadata.getClusterCount() << " clusters " << std::endl;
    fragmentStorage_.resize(tileMetadata.getClusterCount());
    ISAAC_THREAD_CERR << "Resizing fragment storage done for " <<  tileMetadata.getClusterCount() << " clusters " << std::endl;

    const MatchTally::FileTallyList &fileTallyList = matchTally.getFileTallyList(tileMetadata);
    std::vector<Match>::const_iterator barcodeMatchListBegin = matchList.begin();
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList_)
    {
        const unsigned long tileBarcodeMatchCount = std::accumulate(
            fileTallyList.begin(), fileTallyList.end(), 0UL,
            bind(std::plus<unsigned long>(),
                 _1,
                 boost::bind<unsigned long>(&MatchTally::FileTally::getBarcodeMatchCount, _2, barcode.getIndex())));

        if (tileBarcodeMatchCount)
        {
            // we could do determineTemplateLength on multiple threads. The current assumption is that
            // doing it just before using gives some memory cache efficiency benefit which compensates for
            // the absence of parallelization. RP: looks like a wrong assumption though...
            const std::vector<reference::Contig> &barcodeContigList = contigList_.at(barcode.getReferenceIndex());
            const flowcell::ReadMetadataList &tileReads = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex()).getReadMetadataList();
            const RestOfGenomeCorrection restOfGenomeCorrection(barcodeContigList, tileReads);
            TemplateLengthStatistics &templateLengthStatistics = barcodeTemplateLengthStatistics.at(barcode.getIndex());
            if (!templateLengthStatistics.isStable() || perTileTls_)
            {
                ISAAC_THREAD_CERR << "Determining template length for " << tileMetadata << ", " << barcode  << " on " << tileBarcodeMatchCount << " matches." << std::endl;

                templateLengthStatistics =
                    determineTemplateLength(
                        tileMetadata, barcodeContigList, barcodeSequencingAdapters_.at(barcode.getIndex()),
                        barcodeMatchListBegin, barcodeMatchListBegin + tileBarcodeMatchCount,
                        bclData, 0);

                ISAAC_THREAD_CERR << "Determining template length done for " << tileMetadata << ", " << barcode << ":" << templateLengthStatistics << std::endl;
            }
            else
            {
                ISAAC_THREAD_CERR << "Using known template length for " << tileMetadata << ", " << barcode  << " on " << tileBarcodeMatchCount << " matches: " << templateLengthStatistics << std::endl;
            }

            threadStats_[0].recordTemplateLengthStatistics(barcode, templateLengthStatistics);

            ISAAC_THREAD_CERR << "Selecting matches on " <<  computeThreads_.size() << " threads for " <<
                tileMetadata << "," << barcode << std::endl;
            computeThreads_.execute(boost::bind(&MatchSelector::processMatchList, this,
                                                boost::ref(barcodeContigList),
                                                boost::ref(restOfGenomeCorrection),
                                                boost::ref(barcodeSequencingAdapters_.at(barcode.getIndex())),
                                                std::make_pair(barcodeMatchListBegin, barcodeMatchListBegin + tileBarcodeMatchCount),
                                                boost::ref(tileMetadata), boost::ref(bclData),
                                                boost::cref(templateLengthStatistics), _1));

            ISAAC_THREAD_CERR << "Selecting matches done on " <<  computeThreads_.size() << " threads for " <<
                tileMetadata  << "," << barcode << std::endl;
        }

        barcodeMatchListBegin += tileBarcodeMatchCount;
    }
    ISAAC_ASSERT_MSG(matchList.end() == barcodeMatchListBegin, "Expected to reach the end of the tile match list");

    BOOST_FOREACH(const matchSelector::MatchSelectorStats &threadStats, threadStats_)
    {
        allStats_.at(tileMetadata.getIndex()) += threadStats;
    }
}

} // namespace alignemnt
} // namespace isaac
