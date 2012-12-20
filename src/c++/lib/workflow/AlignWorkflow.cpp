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
#include "alignment/matchSelector/BufferingFragmentStorage.hh"
#include "alignment/MatchFinder.hh"
#include "alignment/MatchSelector.hh"
#include "alignment/SeedLoader.hh"
#include "build/Build.hh"
#include "casava/CasavaIntegration.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/FileSystem.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "io/FragmentReader.hh"
#include "reports/AlignmentReportGenerator.hh"

namespace isaac
{
namespace workflow
{

AlignWorkflow::AlignWorkflow(
    const std::vector<std::string> &argv,
    const std::string &casavaArgv,
    const std::vector<flowcell::Layout> &flowcellLayoutList,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const bool allowVariableFastqLength,
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
    const bool putUnalignedInTheBack,
    const bool clipSemialigned,
    const unsigned gappedMismatchesMax,
    const bool scatterRepeats,
    const int gapMatchScore,
    const int gapMismatchScore,
    const int gapOpenScore,
    const int gapExtendScore,
    const alignment::TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore,
    const unsigned inputLoadersMax,
    const unsigned tempSaversMax,
    const unsigned tempLoadersMax,
    const unsigned outputSaversMax,
    const build::GapRealignerMode realignGaps,
    const int bamGzipLevel,
    const double expectedBgzfCompressionRatio,
    const bool keepDuplicates,
    const std::string &binRegexString,
    const common::ScoopedMallocBlock::Mode memoryControl,
    const std::vector<size_t> &clusterIdList,
    const alignment::TemplateLengthStatistics &userTemplateLengthStatistics,
    const reports::AlignmentReportGenerator::ImageFileFormat statsImageFormat)
    : argv_(argv)
    , casavaArgv_(casavaArgv)
    , flowcellLayoutList_(flowcellLayoutList)
    , tempDirectory_(tempDirectory)
    , statsDirectory_(outputDirectory/"Stats")
    , reportsDirectory_(outputDirectory/"Reports")
    , projectsDirectory_(outputDirectory/"Projects")
    , matchSelectorStatsXmlPath_(statsDirectory_ / "MatchSelectorStats.xml")
    , coresMax_(maxThreadCount)
    , repeatThreshold_(repeatThreshold)
    , mateDriftRange_(mateDriftRange)
    , barcodeMetadataList_(barcodeMetadataList)
    , allowVariableFastqLength_(allowVariableFastqLength)
    , ignoreMissingBcls_(ignoreMissingBcls)
    , ignoreMissingFilters_(ignoreMissingFilters)
    , firstPassSeeds_(firstPassSeeds)
    , matchesPerBin_(matchesPerBin)
    , availableMemory_(availableMemory)
    , mapqThreshold_(mapqThreshold)
    , pfOnly_(pfOnly)
    , baseQualityCutoff_(baseQualityCutoff)
    , keepUnaligned_(keepUnaligned)
    , putUnalignedInTheBack_(putUnalignedInTheBack)
    , clipSemialigned_(clipSemialigned)
    , gappedMismatchesMax_(gappedMismatchesMax)
    , scatterRepeats_(scatterRepeats)
    , gapMatchScore_(gapMatchScore)
    , gapMismatchScore_(gapMismatchScore)
    , gapOpenScore_(gapOpenScore)
    , gapExtendScore_(gapExtendScore)
    , dodgyAlignmentScore_(dodgyAlignmentScore)
    , inputLoadersMax_(inputLoadersMax)
    , tempSaversMax_(tempSaversMax)
    , tempLoadersMax_(tempLoadersMax)
    , outputSaversMax_(outputSaversMax)
    , realignGaps_(realignGaps)
    , bamGzipLevel_(bamGzipLevel)
    , expectedBgzfCompressionRatio_(expectedBgzfCompressionRatio)
    , keepDuplicates_(keepDuplicates)
    , binRegexString_(binRegexString)
    , memoryControl_(memoryControl)
    , userTemplateLengthStatistics_(userTemplateLengthStatistics)
    , demultiplexingStatsXmlPath_(statsDirectory_ / "DemultiplexingStats.xml")
    , statsImageFormat_(statsImageFormat)
    , sortedReferenceXmlList_(loadSortedReferenceXml(referenceMetadataList))
    , state_(Start)
    , findMatchesTransition_(flowcellLayoutList_,
                             barcodeMetadataList_,
                             allowVariableFastqLength_,
                             ignoreMissingBcls_,
                             firstPassSeeds_,
                             tempDirectory_,
                             demultiplexingStatsXmlPath_,
                             coresMax_,
                             repeatThreshold_,
                             neighborhoodSizeThreshold,
                             ignoreNeighbors,
                             ignoreRepeats,
                             inputLoadersMax_,
                             tempSaversMax_,
                             memoryControl_,
                             clusterIdList,
                             sortedReferenceXmlList_)
      // dummy initialization. Will be replaced with real object once match finding is over
    , foundMatchesMetadata_(tempDirectory_, barcodeMetadataList_, 0, sortedReferenceXmlList_)
{
    const std::vector<bfs::path> createList = boost::assign::list_of
        (tempDirectory_)(outputDirectory)(statsDirectory_)(reportsDirectory_)(projectsDirectory_);
    common::createDirectories(createList);

    BOOST_FOREACH(const flowcell::Layout &layout, flowcellLayoutList)
    {
        ISAAC_THREAD_CERR << "Aligner: adding base-calls directory " << layout.getBaseCallsDirectory() << std::endl;
    }
}

reference::SortedReferenceXmlList AlignWorkflow::loadSortedReferenceXml(
    const reference::ReferenceMetadataList &referenceMetadataList)
{
    reference::SortedReferenceXmlList ret(referenceMetadataList.size());
    BOOST_FOREACH(const reference::ReferenceMetadata &reference, referenceMetadataList)
    {
        const unsigned referenceIndex = &reference - &referenceMetadataList.front();
        ret.at(referenceIndex) = reference::loadSortedReferenceXml(reference.getXmlPath());
    }
    return ret;
}

alignWorkflow::FoundMatchesMetadata AlignWorkflow::findMatches() const
{
    return findMatchesTransition_.perform();
}

AlignWorkflow::SelectedMatchesMetadata AlignWorkflow::selectMatches() const
{
    // Assume the vast majority of fragments are distributed according
    // to the first pass seed match distribution.
    const unsigned long matchesPerBin = matchesPerBin_
        ? matchesPerBin_
        : build::Build::estimateOptimumFragmentsPerBin(flowcellLayoutList_, availableMemory_, expectedBgzfCompressionRatio_,
                                                       coresMax_) * firstPassSeeds_;


    alignment::matchSelector::BufferingFragmentStorage fragmentStorage(keepUnaligned_, tempSaversMax_, coresMax_,
                       foundMatchesMetadata_.matchDistribution_, matchesPerBin, tempDirectory_,
                       flowcellLayoutList_, barcodeMetadataList_,
                       flowcell::getMaxTileClulsters(foundMatchesMetadata_.tileMetadataList_),
                       foundMatchesMetadata_.tileMetadataList_.size());

    workflow::alignWorkflow::SelectMatchesTransition transition(
                                fragmentStorage,
                                sortedReferenceXmlList_, tempDirectory_, coresMax_,
                                foundMatchesMetadata_.tileMetadataList_, barcodeMetadataList_,
                                flowcellLayoutList_, repeatThreshold_, mateDriftRange_,
                                allowVariableFastqLength_,
                                ignoreMissingBcls_, ignoreMissingFilters_,
                                inputLoadersMax_, tempLoadersMax_, tempSaversMax_,
                                foundMatchesMetadata_.matchTally_,
                                userTemplateLengthStatistics_, mapqThreshold_, pfOnly_, baseQualityCutoff_,
                                keepUnaligned_, clipSemialigned_,
                                gappedMismatchesMax_, scatterRepeats_,
                                gapMatchScore_, gapMismatchScore_, gapOpenScore_, gapExtendScore_, dodgyAlignmentScore_);

    ISAAC_THREAD_CERR << "Selecting matches using " << matchesPerBin << " matches per bin limit" << std::endl;

    transition.selectMatches(memoryControl_, matchSelectorStatsXmlPath_);

    ISAAC_THREAD_CERR << "Selecting matches done using " << matchesPerBin << " matches per bin limit" << std::endl;

    return transition.getBinMetadata();
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
                       sortedReferenceXmlList_,
                       projectsDirectory_,
                       tempLoadersMax_, coresMax_, outputSaversMax_, realignGaps_,
                       bamGzipLevel_, expectedBgzfCompressionRatio_, keepDuplicates_, clipSemialigned_,
                       binRegexString_, alignment::TemplateBuilder::Unknown == dodgyAlignmentScore_,
                       keepUnaligned_, putUnalignedInTheBack_);
    {
        common::ScoopedMallocBlock  mallocBlock(memoryControl_);
        build.run(mallocBlock);
    }
    build.dumpStats(statsDirectory_ / "BuildStats.xml");
    ISAAC_THREAD_CERR << "Generating the BAM files done" << std::endl;
    return build.getBarcodeBamMapping();
}

void AlignWorkflow::resetCasava(const build::BarcodeBamMapping & barcodeBamMapping) const
{
#ifdef HAVE_CASAVA
    ISAAC_THREAD_CERR << "Configuring CASAVA variant calling" << std::endl;
    casava::CasavaIntegration casavaIntegration(tempDirectory_, projectsDirectory_,
                                                barcodeMetadataList_, barcodeBamMapping,
                                                sortedReferenceXmlList_, flowcellLayoutList_);
    casavaIntegration.reset(casavaArgv_);
    ISAAC_THREAD_CERR << "Configuring CASAVA variant calling done" << std::endl;
#endif
}

void AlignWorkflow::resumeCasava(const build::BarcodeBamMapping & barcodeBamMapping) const
{
#ifdef HAVE_CASAVA
    ISAAC_THREAD_CERR << "Executing CASAVA variant calling" << std::endl;
    casava::CasavaIntegration casavaIntegration(tempDirectory_, projectsDirectory_, barcodeMetadataList_,
                                                barcodeBamMapping, sortedReferenceXmlList_, flowcellLayoutList_);
    casavaIntegration.execute(coresMax_);
    ISAAC_THREAD_CERR << "Executing CASAVA variant calling done" << std::endl;
#endif
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
    step();
    ISAAC_ASSERT_MSG(CasavaResetDone == state_, "Unexpected state");
    step();
    ISAAC_ASSERT_MSG(CasavaResumeDone == state_, "Unexpected state");
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
    case BamDone:
    {
        return CasavaResetDone;
    }
    case CasavaResetDone:
    {
        return CasavaResumeDone;
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
        alignWorkflow::FoundMatchesMetadata foundMatchesMetadata = findMatches();
        swap(foundMatchesMetadata_, foundMatchesMetadata);
        state_ = getNextState();
        break;
    }
    case MatchFinderDone:
    {
        SelectedMatchesMetadata selectedMatchesMetadata = selectMatches();
        swap(selectedMatchesMetadata_, selectedMatchesMetadata);
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
    case BamDone:
    {
        resetCasava(barcodeBamMapping_);
        state_ = getNextState();
        break;
    }
    case CasavaResetDone:
    {
        resumeCasava(barcodeBamMapping_);
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

AlignWorkflow::State AlignWorkflow::rewind(AlignWorkflow::State to)
{
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
    case CasavaResetDone:
    {
        if (Start == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from Start to CasavaResetDone is not possible"));}
        if (MatchFinderDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchFinderDone to CasavaResetDone is not possible"));}
        if (MatchSelectorDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchSelectorDone to CasavaResetDone is not possible"));}
        if (AlignmentReportsDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from AlignmentReportsDone to CasavaResetDone is not possible"));}
        if (BamDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from BamDone to CasavaResetDone is not possible"));}
        state_ = CasavaResetDone;
        ISAAC_THREAD_CERR << "Workflow state rewind to CasavaResetDone successful" << std::endl;
        break;
    }
    case CasavaResumeDone:
    {
        if (Start == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from Start to CasavaResumeDone is not possible"));}
        if (MatchFinderDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchFinderDone to CasavaResumeDone is not possible"));}
        if (MatchSelectorDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchSelectorDone to CasavaResumeDone is not possible"));}
        if (AlignmentReportsDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchSelectorDone to CasavaResumeDone is not possible"));}
        if (BamDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from BamDone to CasavaResumeDone is not possible"));}
        if (CasavaResetDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from CasavaResetDone to CasavaResumeDone is not possible"));}
        state_ = CasavaResumeDone;
        ISAAC_THREAD_CERR << "Workflow state rewind to CasavaResumeDone successful" << std::endl;
        break;
    }
    default:
    {
        assert(false);
        break;
    }
    }

    return state_;
}


} // namespace workflow
} // namespace isaac
