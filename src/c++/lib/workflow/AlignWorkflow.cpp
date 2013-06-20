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
 ** \file AlignWorkflow.cpp
 **
 ** \brief see AlignWorkflow.hh
 **
 ** \author Come Raczy
 **/

#include <algorithm>
#include <numeric>
#include <cassert>
#include <cstring>
#include <cerrno>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>

#include "workflow/AlignWorkflow.hh"
#include "alignment/matchSelector/BinningFragmentStorage.hh"
#include "alignment/MatchFinder.hh"
#include "alignment/MatchSelector.hh"
#include "alignment/SeedLoader.hh"
#include "build/Build.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/FileSystem.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "reports/AlignmentReportGenerator.hh"

namespace isaac
{
namespace workflow
{

AlignWorkflow::AlignWorkflow(
    const std::vector<std::string> &argv,
    const std::vector<flowcell::Layout> &flowcellLayoutList,
    const unsigned seedLength,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const bool allowVariableFastqLength,
    const bool cleanupIntermediary,
    const bool ignoreMissingBcls,
    const bool ignoreMissingFilters,
    const unsigned firstPassSeeds,
    const unsigned long matchesPerBin,
    const reference::ReferenceMetadataList &referenceMetadataList,
    const bfs::path &tempDirectory,
    const bfs::path &outputDirectory,
    const unsigned int maxThreadCount,
    const unsigned repeatThreshold,
    const int mateDriftRange,
    const unsigned neighborhoodSizeThreshold,
    const unsigned long availableMemory,
    const bool ignoreNeighbors,
    const bool ignoreRepeats,
    const unsigned mapqThreshold,
    const bool pfOnly,
    const unsigned baseQualityCutoff,
    const bool keepUnaligned,
    const bool preSortBins,
    const bool putUnalignedInTheBack,
    const bool realignGapsVigorously,
    const bool realignDodgyFragments,
    const unsigned realignedGapsPerFragment,
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
    const unsigned inputLoadersMax,
    const unsigned tempSaversMax,
    const unsigned tempLoadersMax,
    const unsigned outputSaversMax,
    const build::GapRealignerMode realignGaps,
    const int bamGzipLevel,
    const std::vector<std::string> &bamHeaderTags,
    const double expectedBgzfCompressionRatio,
    const bool singleLibrarySamples,
    const bool keepDuplicates,
    const bool markDuplicates,
    const std::string &binRegexString,
    const common::ScoopedMallocBlock::Mode memoryControl,
    const std::vector<std::size_t> &clusterIdList,
    const alignment::TemplateLengthStatistics &userTemplateLengthStatistics,
    const reports::AlignmentReportGenerator::ImageFileFormat statsImageFormat,
    const bool bufferBins,
    const bool qScoreBin,
    const boost::array<char, 256> &fullBclQScoreTable,
    const OptionalFeatures optionalFeatures)
    : argv_(argv)
    , flowcellLayoutList_(flowcellLayoutList)
    , seedLength_(seedLength)
    , tempDirectory_(tempDirectory)
    , statsDirectory_(outputDirectory/"Stats")
    , reportsDirectory_(outputDirectory/"Reports")
    , projectsDirectory_(outputDirectory/"Projects")
    , matchSelectorStatsXmlPath_(statsDirectory_ / "MatchSelectorStats.xml")
    , coresMax_(maxThreadCount)
    , repeatThreshold_(repeatThreshold)
    , mateDriftRange_(mateDriftRange)
    , neighborhoodSizeThreshold_(neighborhoodSizeThreshold)
    , ignoreNeighbors_(ignoreNeighbors)
    , ignoreRepeats_(ignoreRepeats)
    , clusterIdList_(clusterIdList)
    , barcodeMetadataList_(barcodeMetadataList)
    , allowVariableFastqLength_(allowVariableFastqLength)
    , cleanupIntermediary_(cleanupIntermediary)
    , ignoreMissingBcls_(ignoreMissingBcls)
    , ignoreMissingFilters_(ignoreMissingFilters)
    , firstPassSeeds_(firstPassSeeds)
    , matchesPerBin_(matchesPerBin)
    , availableMemory_(availableMemory)
    , mapqThreshold_(mapqThreshold)
    , pfOnly_(pfOnly)
    , baseQualityCutoff_(baseQualityCutoff)
    , keepUnaligned_(keepUnaligned)
    , preSortBins_(preSortBins)
    , putUnalignedInTheBack_(putUnalignedInTheBack)
    , realignGapsVigorously_(realignGapsVigorously)
    , realignDodgyFragments_(realignDodgyFragments)
    , realignedGapsPerFragment_(realignedGapsPerFragment)
    , clipSemialigned_(clipSemialigned)
    , clipOverlapping_(clipOverlapping)
    , scatterRepeats_(scatterRepeats)
    , gappedMismatchesMax_(gappedMismatchesMax)
    , avoidSmithWaterman_(avoidSmithWaterman)
    , gapMatchScore_(gapMatchScore)
    , gapMismatchScore_(gapMismatchScore)
    , gapOpenScore_(gapOpenScore)
    , gapExtendScore_(gapExtendScore)
    , minGapExtendScore_(minGapExtendScore)
    , semialignedGapLimit_(semialignedGapLimit)
    , dodgyAlignmentScore_(dodgyAlignmentScore)
    , inputLoadersMax_(inputLoadersMax)
    , tempSaversMax_(tempSaversMax)
    , tempLoadersMax_(tempLoadersMax)
    , outputSaversMax_(outputSaversMax)
    , realignGaps_(realignGaps)
    , bamGzipLevel_(bamGzipLevel)
    , bamHeaderTags_(bamHeaderTags)
    , expectedBgzfCompressionRatio_(expectedBgzfCompressionRatio)
    , singleLibrarySamples_(singleLibrarySamples)
    , keepDuplicates_(keepDuplicates)
    , markDuplicates_(markDuplicates)
    , bufferBins_(bufferBins)
    , qScoreBin_(qScoreBin)
    , fullBclQScoreTable_(fullBclQScoreTable)
    , optionalFeatures_(optionalFeatures)
    , binRegexString_(binRegexString)
    , memoryControl_(memoryControl)
    , userTemplateLengthStatistics_(userTemplateLengthStatistics)
    , demultiplexingStatsXmlPath_(statsDirectory_ / "DemultiplexingStats.xml")
    , statsImageFormat_(statsImageFormat)
    , sortedReferenceMetadataList_(loadSortedReferenceXml(seedLength, referenceMetadataList))
    , state_(Start)
      // dummy initialization. Will be replaced with real object once match finding is over
    , foundMatchesMetadata_(tempDirectory_, barcodeMetadataList_, 0, sortedReferenceMetadataList_)
{
    const std::vector<bfs::path> createList = boost::assign::list_of
        (tempDirectory_)(outputDirectory)(statsDirectory_)(reportsDirectory_)(projectsDirectory_);
    common::createDirectories(createList);

    BOOST_FOREACH(const flowcell::Layout &layout, flowcellLayoutList)
    {
        ISAAC_THREAD_CERR << "Aligner: adding base-calls path " << layout.getBaseCallsPath() << std::endl;
    }
}

reference::SortedReferenceMetadataList AlignWorkflow::loadSortedReferenceXml(
    const unsigned seedLength,
    const reference::ReferenceMetadataList &referenceMetadataList)
{
    reference::SortedReferenceMetadataList ret(referenceMetadataList.size());
    BOOST_FOREACH(const reference::ReferenceMetadata &reference, referenceMetadataList)
    {
        const unsigned referenceIndex = &reference - &referenceMetadataList.front();
        reference::SortedReferenceMetadata &ref = ret.at(referenceIndex);
        ref = reference::loadSortedReferenceXml(reference.getXmlPath());
        if (!ref.supportsSeedLength(seedLength))
        {
            boost::format msg = boost::format("Sorted reference %s does not support seed length %d") %
                reference.getXmlPath().string() % seedLength;
            BOOST_THROW_EXCEPTION(common::PreConditionException(msg.str()));
        }
    }
    return ret;
}

void AlignWorkflow::findMatches(alignWorkflow::FoundMatchesMetadata &foundMatches) const
{
    alignWorkflow::FindMatchesTransition findMatchesTransition(
        flowcellLayoutList_,
        barcodeMetadataList_,
        allowVariableFastqLength_,
        ignoreMissingBcls_,
        firstPassSeeds_,
        availableMemory_,
        tempDirectory_,
        demultiplexingStatsXmlPath_,
        coresMax_,
        repeatThreshold_,
        neighborhoodSizeThreshold_,
        ignoreNeighbors_,
        ignoreRepeats_,
        inputLoadersMax_,
        tempSaversMax_,
        memoryControl_,
        clusterIdList_,
        sortedReferenceMetadataList_);

    if (16 == seedLength_)
    {
        findMatchesTransition.perform<isaac::oligo::ShortKmerType>(foundMatches);
    }
    else if (32 == seedLength_)
    {
        findMatchesTransition.perform<oligo::KmerType>(foundMatches);
    }
    else if (64 == seedLength_)
    {
        findMatchesTransition.perform<oligo::LongKmerType>(foundMatches);
    }
    else
    {
        ISAAC_ASSERT_MSG(false, "Unexpected seed length " << seedLength_);
    }
}

void AlignWorkflow::cleanupMatches() const
{
    ISAAC_THREAD_CERR << "Removing intermediary match files" << std::endl;
    unsigned removed = 0;
    BOOST_FOREACH(const flowcell::TileMetadata &tile,  foundMatchesMetadata_.tileMetadataList_)
    {
        BOOST_FOREACH(const alignment::MatchTally::FileTally &file, foundMatchesMetadata_.matchTally_.getFileTallyList(tile))
        {
            removed += boost::filesystem::remove(file.path_);
        }
    }
    ISAAC_THREAD_CERR << "Removing intermediary match files done. " << removed << " files removed." << std::endl;
}


void AlignWorkflow::selectMatches(
    alignment::matchSelector::FragmentStorage &fragmentStorage) const
{
    ISAAC_TRACE_STAT("AlignWorkflow::selectMatches fragmentStorage ")

    workflow::alignWorkflow::SelectMatchesTransition transition(
        (bufferBins_ ? 3 : 2),
        fragmentStorage, foundMatchesMetadata_.matchDistribution_,
        sortedReferenceMetadataList_, tempDirectory_, coresMax_,
        foundMatchesMetadata_.tileMetadataList_, barcodeMetadataList_,
        flowcellLayoutList_, repeatThreshold_, mateDriftRange_,
        allowVariableFastqLength_,
        ignoreMissingBcls_, ignoreMissingFilters_,
        inputLoadersMax_, tempLoadersMax_, tempSaversMax_,
        foundMatchesMetadata_.matchTally_,
        userTemplateLengthStatistics_, mapqThreshold_, pfOnly_, baseQualityCutoff_,
        keepUnaligned_, clipSemialigned_, clipOverlapping_,
        scatterRepeats_, gappedMismatchesMax_, avoidSmithWaterman_,
        gapMatchScore_, gapMismatchScore_, gapOpenScore_, gapExtendScore_, minGapExtendScore_, semialignedGapLimit_,
        dodgyAlignmentScore_, qScoreBin_, fullBclQScoreTable_, optionalFeatures_ & BamZX);

    transition.selectMatches(memoryControl_, matchSelectorStatsXmlPath_);
}

void AlignWorkflow::cleanupBins() const
{
    ISAAC_THREAD_CERR << "Removing intermediary bin files" << std::endl;
    unsigned removed = 0;
    BOOST_FOREACH(const alignment::BinMetadata &bin, selectedMatchesMetadata_)
    {
        removed += boost::filesystem::remove(bin.getPath());
    }
    ISAAC_THREAD_CERR << "Removing intermediary bin files done. " << removed << " files removed." << std::endl;
}

void AlignWorkflow::selectMatches(SelectedMatchesMetadata &binPaths) const
{
    // Assume the vast majority of fragments are distributed according
    // to the first pass seed match distribution.
    const unsigned long matchesPerBin = matchesPerBin_
        ? matchesPerBin_
        : build::Build::estimateOptimumFragmentsPerBin(flowcellLayoutList_, availableMemory_, expectedBgzfCompressionRatio_,
                                                       coresMax_) * firstPassSeeds_;

    ISAAC_TRACE_STAT("AlignWorkflow::selectMatches ")

    if (!bufferBins_)
    {
        alignment::matchSelector::BinningFragmentStorage fragmentStorage(
            keepUnaligned_, preSortBins_, tempSaversMax_, coresMax_,
            foundMatchesMetadata_.matchDistribution_, matchesPerBin, tempDirectory_,
            flowcellLayoutList_, barcodeMetadataList_,
            flowcell::getMaxTileClusters(foundMatchesMetadata_.tileMetadataList_),
            foundMatchesMetadata_.tileMetadataList_.size());

        ISAAC_THREAD_CERR << "Selecting matches using " << matchesPerBin << " matches per bin limit" << std::endl;
        selectMatches(fragmentStorage);
        AlignWorkflow::SelectedMatchesMetadata ret;
        fragmentStorage.close(binPaths);
    }
    else
    {
        alignment::matchSelector::BufferingFragmentStorage fragmentStorage(
            keepUnaligned_, preSortBins_, tempSaversMax_, coresMax_,
            foundMatchesMetadata_.matchDistribution_, matchesPerBin, tempDirectory_,
            flowcellLayoutList_, barcodeMetadataList_,
            flowcell::getMaxTileClusters(foundMatchesMetadata_.tileMetadataList_),
            foundMatchesMetadata_.tileMetadataList_.size(),
            "skip-empty" == binRegexString_);

        ISAAC_THREAD_CERR << "Selecting matches using " << matchesPerBin << " matches per bin limit" << std::endl;
        selectMatches(fragmentStorage);
        AlignWorkflow::SelectedMatchesMetadata ret;
        fragmentStorage.close(binPaths);
    }
    ISAAC_THREAD_CERR << "Selecting matches done using " << matchesPerBin << " matches per bin limit. Produced " << binPaths.size() << " bins." << std::endl;
}

void AlignWorkflow::generateAlignmentReports() const
{
    ISAAC_THREAD_CERR << "Generating the match selector reports from " << matchSelectorStatsXmlPath_ << std::endl;
    reports::AlignmentReportGenerator reportGenerator(flowcellLayoutList_, barcodeMetadataList_,
                                                  matchSelectorStatsXmlPath_, demultiplexingStatsXmlPath_,
                                                  tempDirectory_, reportsDirectory_,
                                                  statsImageFormat_);
    reportGenerator.run();
    ISAAC_THREAD_CERR << "Generating the match selector reports done from " << matchSelectorStatsXmlPath_ << std::endl;
}

const build::BarcodeBamMapping AlignWorkflow::generateBam(const SelectedMatchesMetadata &binPaths) const
{
    ISAAC_THREAD_CERR << "Generating the BAM files" << std::endl;

    build::Build build(argv_, flowcellLayoutList_, foundMatchesMetadata_.tileMetadataList_, barcodeMetadataList_,
                       binPaths,
                       sortedReferenceMetadataList_,
                       projectsDirectory_,
                       tempLoadersMax_, coresMax_, outputSaversMax_, realignGaps_,
                       bamGzipLevel_, bamHeaderTags_, expectedBgzfCompressionRatio_, singleLibrarySamples_,
                       keepDuplicates_, markDuplicates_,
                       realignGapsVigorously_, realignDodgyFragments_, realignedGapsPerFragment_,
                       clipSemialigned_, binRegexString_,
                       alignment::TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED == dodgyAlignmentScore_ ?
                           0 : boost::numeric_cast<unsigned char>(dodgyAlignmentScore_),
                       keepUnaligned_, putUnalignedInTheBack_,
                       build::IncludeTags(
                           optionalFeatures_ & BamAS,
                           optionalFeatures_ & BamBC,
                           optionalFeatures_ & BamNM,
                           optionalFeatures_ & BamOC,
                           optionalFeatures_ & BamRG,
                           optionalFeatures_ & BamSM,
                           optionalFeatures_ & BamZX,
                           optionalFeatures_ & BamZY));
    {
        common::ScoopedMallocBlock  mallocBlock(memoryControl_);
        build.run(mallocBlock);
    }
    build.dumpStats(statsDirectory_ / "BuildStats.xml");
    ISAAC_THREAD_CERR << "Generating the BAM files done" << std::endl;
    return build.getBarcodeBamMapping();
}

void AlignWorkflow::run()
{
    ISAAC_ASSERT_MSG(Start == state_, "Unexpected state");
    step();
    ISAAC_ASSERT_MSG(MatchFinderDone == state_, "Unexpected state");
    step();
    ISAAC_ASSERT_MSG(MatchSelectorDone == state_, "Unexpected state");
    step();
    ISAAC_ASSERT_MSG(AlignmentReportsDone == state_, "Unexpected state");
    step();
    ISAAC_ASSERT_MSG(BamDone == state_, "Unexpected state");
}

AlignWorkflow::State AlignWorkflow::getNextState() const
{
    switch (state_)
    {
    case Start:
    {
        return MatchFinderDone;
    }
    case MatchFinderDone:
    {
        return MatchSelectorDone;
    }
    case MatchSelectorDone:
    {
        return AlignmentReportsDone;
    }
    case AlignmentReportsDone:
    {
        return BamDone;
    }
    case Finish:
    {
        return Finish;
    }
    default:
    {
        ISAAC_ASSERT_MSG(false, "Invalid state value");
        return Invalid;
    }
    }
}

AlignWorkflow::State AlignWorkflow::step()
{
    using std::swap;
    switch (state_)
    {
    case Start:
    {
        findMatches(foundMatchesMetadata_);
        state_ = getNextState();
        break;
    }
    case MatchFinderDone:
    {
        selectMatches(selectedMatchesMetadata_);
        state_ = getNextState();
        break;
    }
    case MatchSelectorDone:
    {
        generateAlignmentReports();
        state_ = getNextState();
        break;
    }
    case AlignmentReportsDone:
    {
        barcodeBamMapping_ = generateBam(selectedMatchesMetadata_);
        state_ = getNextState();
        break;
    }
    case Finish:
    {
        ISAAC_THREAD_CERR << "Already at the Finish state" << std::endl;
        break;
    }
    default:
    {
        ISAAC_ASSERT_MSG(false, "Invalid state");
        break;
    }
    }
    return state_;
}

void AlignWorkflow::cleanupIntermediary()
{
    switch (state_)
    {
    case Finish:
    {
        cleanupBins();
        //fall through
    }
    case AlignmentReportsDone:
    case MatchSelectorDone:
    {
        cleanupMatches();
        //fall through
    }
    case MatchFinderDone:
    case Start:
    {
        break;
    }
    default:
    {
        ISAAC_ASSERT_MSG(false, "Invalid state " << state_);
        break;
    }
    }
}

/**
 * \return the initial state from which the rewind occurred
 */
AlignWorkflow::State AlignWorkflow::rewind(AlignWorkflow::State to)
{
    AlignWorkflow::State ret = state_;
    switch (to)
    {
    case Last:
    {
        // Nothing to do. We're at the Last state by definition
        break;
    }
    case Start:
    {
        // Start is always possible
        state_ = Start;
        break;
    }
    case MatchFinderDone:
    {
        if (Start == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from Start to MatchFinderDone is not possible"));}
        state_ = MatchFinderDone;
        ISAAC_THREAD_CERR << "Workflow state rewind to MatchFinderDone successful" << std::endl;
        break;
    }
    case MatchSelectorDone:
    {
        if (Start == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from Start to MatchSelectorDone is not possible"));}
        if (MatchFinderDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchFinderDone to MatchSelectorDone is not possible"));}
        state_ = MatchSelectorDone;
        ISAAC_THREAD_CERR << "Workflow state rewind to MatchSelectorDone successful" << std::endl;
        break;
    }
    case AlignmentReportsDone:
    {
        if (Start == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from Start to MatchSelectorReportsDone is not possible"));}
        if (MatchFinderDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchFinderDone to MatchSelectorReportsDone is not possible"));}
        if (MatchSelectorDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchSelectorDone to MatchSelectorReportsDone is not possible"));}
        state_ = AlignmentReportsDone;
        ISAAC_THREAD_CERR << "Workflow state rewind to AlignmentReportsDone successful" << std::endl;
        break;
    }
    case BamDone:
    {
        if (Start == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from Start to BamDone is not possible"));}
        if (MatchFinderDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchFinderDone to BamDone is not possible"));}
        if (MatchSelectorDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchSelectorDone to BamDone is not possible"));}
        if (AlignmentReportsDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchSelectorDone to BamDone is not possible"));}
        state_ = BamDone;
        ISAAC_THREAD_CERR << "Workflow state rewind to BamDone successful" << std::endl;
        break;
    }
    default:
    {
        assert(false);
        break;
    }
    }

    return ret;
}


} // namespace workflow
} // namespace isaac
