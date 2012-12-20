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
 ** \file FindMatchesTransition.hh
 **
 ** \brief Top level component to control the analysis process.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FIND_MATCHES_TRANSITION_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_FIND_MATCHES_TRANSITION_HH

#include "alignment/MatchTally.hh"
#include "alignment/MatchDistribution.hh"
#include "alignment/SeedMetadata.hh"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "common/Threads.hpp"
#include "demultiplexing/BarcodeLoader.hh"
#include "demultiplexing/BarcodeResolver.hh"
#include "demultiplexing/DemultiplexingStats.hh"
#include "flowcell/Layout.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "oligo/Kmer.hh"
#include "oligo/Permutations.hh"
#include "reference/ReferenceMetadata.hh"
#include "reference/SortedReferenceXml.hh"

#include "workflow/alignWorkflow/DataSource.hh"
#include "workflow/alignWorkflow/FoundMatchesMetadata.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

namespace bfs = boost::filesystem;

class FindMatchesTransition: boost::noncopyable
{
public:
    typedef flowcell::TileMetadata TileMetadata;
    typedef flowcell::ReadMetadata ReadMetadata;
    typedef flowcell::ReadMetadataList ReadMetadataList;

    FindMatchesTransition(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const bool allowVariableFastqLength,
        const bool ignoreMissingBcls,
        const unsigned firstPassSeeds,
        const bfs::path &tempDirectory,
        const bfs::path &demultiplexingStatsXmlPath,
        const unsigned maxThreadCount,
        const unsigned repeatThreshold,
        const unsigned neighborhoodSizeThreshold,
        const bool ignoreNeighbors,
        const bool ignoreRepeats,
        const unsigned inputLoadersMax,
        const unsigned tempSaversMax,
        const common::ScoopedMallocBlock::Mode memoryControl,
        const std::vector<size_t> &clusterIdList,
        const reference::SortedReferenceXmlList &sortedReferenceXmlList);

    FoundMatchesMetadata perform() const;

private:
    template<class Archive> friend void serialize(Archive & ar, FindMatchesTransition &, const unsigned int file_version);

    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    const bfs::path tempDirectory_;
    const bfs::path demultiplexingStatsXmlPath_;
    const unsigned coresMax_;
    const unsigned repeatThreshold_;
    const unsigned neighborhoodSizeThreshold_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const bool allowVariableFastqLength_;
    const bool ignoreMissingBcls_;
    const unsigned firstPassSeeds_;
    const bool ignoreNeighbors_;
    const bool ignoreRepeats_;
    const unsigned inputLoadersMax_;
    const unsigned tempSaversMax_;
    const common::ScoopedMallocBlock::Mode memoryControl_;
    const std::vector<size_t> clusterIdList_;

    const reference::SortedReferenceXmlList sortedReferenceXmlList_;

    static const unsigned maxIterations_ = 2;

    static reference::SortedReferenceXmlList loadSortedReferenceXml(
        const reference::ReferenceMetadataList &referenceMetadataList);

    /**
     ** \brief Return the list of seed indexes to use for each iteration.
     **
     ** The first seed of each read is always used for the first iteration. All
     ** subsequent seeds, if any, are used for the second iteration. The outer
     ** vector is for the successive iterations. The inner vector is the list of
     ** seed indexes (in the seedMetadataList_) to use for each iteration.
     **
     ** TODO: implement a more flexible multi-seed management
     **/
    std::vector<std::vector<unsigned> > getSeedIndexListPerIteration(const flowcell::Layout &flowcell) const;

    void resolveBarcodes(
        const flowcell::Layout &flowcell,
        const flowcell::BarcodeMetadataList &barcodeGroup,
        const flowcell::TileMetadataList &allTiles,
        flowcell::TileMetadataList unprocessedTiles,
        alignment::matchFinder::TileClusterInfo &tileClusterInfo,
        demultiplexing::DemultiplexingStats &demultiplexingStats) const;

    flowcell::TileMetadataList findSingleSeedMatches(
        const flowcell::Layout &flowcell,
        const std::vector<unsigned> &seedIndexList,
        flowcell::TileMetadataList &unprocessedTiles,
        alignment::matchFinder::TileClusterInfo &tileClusterInfo,
        DataSource &dataSource,
        demultiplexing::DemultiplexingStats &demultiplexingStats,
        FoundMatchesMetadata &foundMatches) const;

    void findMultiSeedMatches(
        const flowcell::Layout &flowcell,
        const std::vector<unsigned> &seedIndexList,
        flowcell::TileMetadataList &unprocessedTiles,
        alignment::matchFinder::TileClusterInfo &tileClusterInfo,
        DataSource &dataSource,
        FoundMatchesMetadata &foundMatches) const;

    void findLaneMatches(
        const flowcell::Layout &flowcell,
        const unsigned lane,
        const flowcell::BarcodeMetadataList &barcodeGroup,
        flowcell::TileMetadataList &unprocessedTiles,
        DataSource &dataSource,
        demultiplexing::DemultiplexingStats &demultiplexingStats,
        FoundMatchesMetadata &foundMatches) const;

    void processFlowcellTiles(
        const flowcell::Layout& flowcell,
        DataSource &dataSource,
        demultiplexing::DemultiplexingStats &demultiplexingStats,
        FoundMatchesMetadata &foundMatches) const;

    void dumpStats(
        const demultiplexing::DemultiplexingStats &demultiplexingStats,
        const flowcell::TileMetadataList &tileMetadataList) const;

};

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FIND_MATCHES_TRANSITION_HH
