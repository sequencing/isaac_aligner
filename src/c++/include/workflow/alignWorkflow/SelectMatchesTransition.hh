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
 ** \file SelectMatchesTransition.hh
 **
 ** \brief Threading, memory and and file management for converting matches into aligned clusters.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_SELECT_MATCHES_TRANSITION_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_SELECT_MATCHES_TRANSITION_HH

#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "flowcell/ReadMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "alignment/Cluster.hh"
#include "alignment/FragmentBuilder.hh"
#include "alignment/TemplateBuilder.hh"
#include "alignment/Match.hh"
#include "alignment/MatchDistribution.hh"
#include "alignment/MatchTally.hh"
#include "alignment/SeedMetadata.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "alignment/MatchSelector.hh"
#include "alignment/matchSelector/FragmentStorage.hh"
#include "alignment/matchSelector/MatchSelectorStats.hh"
#include "alignment/matchSelector/ParallelMatchLoader.hh"
#include "alignment/matchSelector/SemialignedEndsClipper.hh"
#include "common/Threads.hpp"
#include "flowcell/BarcodeMetadata.hh"
#include "io/FastqLoader.hh"
#include "reference/Contig.hh"
#include "reference/SortedReferenceMetadata.hh"
#include "rta/BclReader.hh"
#include "workflow/alignWorkflow/BamDataSource.hh"
#include "workflow/alignWorkflow/BclBgzfDataSource.hh"
#include "workflow/alignWorkflow/BclDataSource.hh"
#include "workflow/alignWorkflow/FastqDataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

class SelectMatchesTransition: boost::noncopyable
{
public:
    typedef flowcell::TileMetadataList TileMetadataList;
    typedef flowcell::ReadMetadataList ReadMetadataList;
    /// Construction of an instance for a given reference
    SelectMatchesTransition(
        const unsigned ioOverlapParallelization,
        alignment::matchSelector::FragmentStorage &fragmentStorage,
        const alignment::MatchDistribution &matchDistribution,
        const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const boost::filesystem::path &tempDirectory,
        unsigned int maxThreadCount,
        const TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned repeatThreshold,
        const int mateDriftRange,
        const bool allowVariableFastqLength,
        const bool ignoreMissingBcls,
        const bool ignoreMissingFilters,
        const unsigned inputLoadersMax,
        const unsigned tempLoadersMax,
        const unsigned tempSaversMax,
        const alignment::MatchTally &matchTally,
        const alignment::TemplateLengthStatistics &defaultTemplateLengthStatistics,
        const unsigned mapqThreshold,
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
        const alignment::TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore,
        const bool qScoreBin,
        const boost::array<char, 256> &fullBclQScoreTable,
        const bool extractClusterXy);

    /**
     ** \brief Select the best match for each of the cluster in the given tile
     ** feed them to the fragmentDispatcher.
     **/
    void selectMatches(
        const common::ScoopedMallocBlock::Mode memoryControl,
        const boost::filesystem::path &matchSelectorStatsXmlPath,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics);

private:
    common::ThreadVector matchLoadThreads_;
    common::ThreadVector inputLoaderThreads_;

    const TileMetadataList tileMetadataList_;
    /**
     * \brief threadBclFilePaths_ gets resized for every tile total readlength. If the tile read lengths
     *        changes from lower to bigger, more threadBclFilePaths_ strings get allocated which breaks the whole
     *        concept of allocating things once. For now this list contains tiles in the processing order so
     *        that the total read length goes only down. TODO: cleanup this mess for example by creating
     *        SelectMatchesTransition only for the group of tiles that have the same geometry.
     */
    const TileMetadataList processOrderTileMetadataList_;
    const flowcell::FlowcellLayoutList flowcellLayoutList_;

    common::ThreadVector ioOverlapThreads_;
    TileMetadataList::const_iterator nextUnprocessedTile_;


    mutable bool loadSlotAvailable_;
    mutable bool flushSlotAvailable_;
    mutable bool computeSlotAvailable_;
    mutable boost::condition_variable stateChangedCondition_;
    mutable boost::mutex slotMutex_;

    const std::vector<alignment::matchSelector::SequencingAdapterList> barcodeSequencingAdapters_;

    const alignment::MatchTally &matchTally_;
    std::vector<std::vector<alignment::Match> > threadMatches_;

    alignment::matchSelector::FragmentStorage &fragmentStorage_;

    alignment::matchSelector::ParallelMatchLoader matchLoader_;
    std::vector<alignment::BclClusters> threadBclData_;
    boost::scoped_ptr<BclBaseCallsSource> bclBaseCallsSource_;
    boost::scoped_ptr<FastqBaseCallsSource> fastqBaseCallsSource_;
    boost::scoped_ptr<BamBaseCallsSource> bamBaseCallsSource_;
    boost::scoped_ptr<BclBgzfBaseCallsSource> bclBgzfBaseCallsSource_;

    alignment::MatchSelector matchSelector_;
    bool qScoreBin_;
    const boost::array<char, 256> &fullBclQScoreTable_;
    bool forceTermination_;

    void acquireSlot(bool &slotAvailable) const
    {
        boost::unique_lock<boost::mutex> lock(slotMutex_);
        while(!slotAvailable){stateChangedCondition_.wait(lock);}
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }
        slotAvailable = false;
    }

    void releaseSlot(bool &slotAvailable, const bool exceptionUnwinding)
    {
        boost::unique_lock<boost::mutex> lock(slotMutex_);
        assert(!slotAvailable && "The slot is must be acquired for it to be released");
        slotAvailable = true;
        if (exceptionUnwinding)
        {
            forceTermination_ = true;
        }

        stateChangedCondition_.notify_all(); //let other thread(s) know load slot is available again
    }

    void acquireLoadSlot() const {acquireSlot(loadSlotAvailable_);}
    void acquireComputeSlot() const {acquireSlot(computeSlotAvailable_);}
    void acquireFlushSlot() const {acquireSlot(flushSlotAvailable_);}

    void releaseLoadSlot(const bool exceptionUnwinding) {releaseSlot(loadSlotAvailable_, exceptionUnwinding);}
    void releaseComputeSlot(const bool exceptionUnwinding) {releaseSlot(computeSlotAvailable_, exceptionUnwinding);}
    void releaseFlushSlot(const bool exceptionUnwinding) {releaseSlot(flushSlotAvailable_, exceptionUnwinding);}

    void processMatchList(
        const std::vector<reference::Contig> &barcodeContigList,
        const alignment::matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<alignment::Match>::const_iterator ourMatchListBegin,
        const std::vector<alignment::Match>::const_iterator ourMatchListEnd,
        const flowcell::TileMetadata & tileMetadata,
        const alignment::BclClusters &bclData,
        const alignment::TemplateLengthStatistics & templateLengthStatistics,
        const unsigned threadNumber);

    void parallelSelect(
        const flowcell::TileMetadata &tileMetadata,
        std::vector<alignment::Match> &matchList,
        const alignment::BclClusters &bclData);

    void loadClusters(const unsigned threadNumber, const flowcell::TileMetadata &tileMetadata);

    /**
     * \brief Processes a single tile
     **/
    void selectTileMatches(
        const unsigned threadNumber,
        const alignment::MatchTally &matchTally,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        common::ScoopedMallocBlock &mallocBlock);

    /**
     ** \brief load all the data from the given tile into the selected destination
     **
     ** The layout and format of the data is:
     **   - 1 byte/base, as in the bcl files (bits are qqqqqqbb)
     **   - bases in a read are contiguous in memory
     **   - reads from a single cluster are contiguous
     **
     ** Note: loading is done for a complete time at a time, so that
     **   - loading a tile can be cached while processing the previous tile
     **   - high level of parallelization won't stumble on open file limit
     **
     ** Note: RVO should take care of the spurious copy. Otherwise, use
     ** boost::move or swap something like that to allow the creation of a const
     **/
    void bclToClusters(
        const rta::ParallelBclMapper<rta::BclReader> &bclMapper,
        const io::FiltersMapper &filtersMapper,
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData) const;

    /**
     * \brief Construct the contig list from the SortedReference XML
     */
    std::vector<reference::Contig> getContigList(
        const reference::SortedReferenceMetadata &sortedReferenceMetadata) const;

    unsigned long getMaxTileMatches(const alignment::MatchTally &matchTally) const
    {
        unsigned long ret = 0;
        BOOST_FOREACH(const flowcell::TileMetadata &tileMetadata, tileMetadataList_)
        {
            const std::vector<alignment::MatchTally::FileTally> &tileFileTally = matchTally.getFileTallyList(tileMetadata);
            ret = std::max(ret,
                           std::accumulate(tileFileTally.begin(), tileFileTally.end(), 0ul,
                                   boost::bind(std::plus<unsigned long>(), _1,
                                               boost::bind(&alignment::MatchTally::FileTally::matchCount_, _2))));
        }
        return ret;
    }

    /**
     * \return Iterator to the first match of the next cluster
     */
    static std::vector<alignment::Match>::const_iterator findNextCluster(
        std::vector<alignment::Match>::const_iterator currentClusterIt,
        std::vector<alignment::Match>::const_iterator endIt)
    {
        if (currentClusterIt == endIt)
        {
            return endIt;
        }
        const unsigned long currentClusterId = currentClusterIt->getCluster();
        const unsigned long currentTileBarcode = currentClusterIt->getTileBarcode();
        ISAAC_THREAD_CERR_DEV_TRACE((boost::format("    match: %s") % *currentClusterIt).str());
        while ((++currentClusterIt != endIt) && (currentClusterId == currentClusterIt->getCluster()))
        {
            ISAAC_THREAD_CERR_DEV_TRACE((boost::format("    match: %s") % *currentClusterIt).str());
            ISAAC_ASSERT_MSG(currentTileBarcode == currentClusterIt->getTileBarcode(), "Matches of the same cluster expected to have the same barcode and tile.");
        }
        return currentClusterIt;
    }
};

} // namespace alignWorkflow
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_SELECT_MATCHES_TRANSITION_HH
