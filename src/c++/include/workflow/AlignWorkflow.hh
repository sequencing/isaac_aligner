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
 ** \file AlignWorkflow.hh
 **
 ** \brief Top level component to controll the analysis process.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_HH

#include <string>
#include <vector>

#include <boost/serialization/access.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>

#include "alignment/BinMetadata.hh"
#include "alignment/MatchTally.hh"
#include "alignment/MatchDistribution.hh"
#include "alignment/SeedMetadata.hh"
#include "alignment/TemplateBuilder.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "build/BarcodeBamMapping.hh"
#include "build/BinSorter.hh"
#include "common/Threads.hpp"
#include "demultiplexing/BarcodeLoader.hh"
#include "demultiplexing/BarcodeResolver.hh"
#include "flowcell/Layout.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "oligo/Kmer.hh"
#include "oligo/Permutations.hh"
#include "reference/ReferenceMetadata.hh"
#include "reference/SortedReferenceXml.hh"

#include "workflow/alignWorkflow/FindMatchesTransition.hh"
#include "workflow/alignWorkflow/SelectMatchesTransition.hh"
#include "workflow/alignWorkflow/FoundMatchesMetadata.hh"

#include "reports/AlignmentReportGenerator.hh"

namespace isaac
{
namespace workflow
{

namespace bfs = boost::filesystem;

class AlignWorkflow: boost::noncopyable
{
public:

    AlignWorkflow(
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
        const unsigned maxThreadCount,
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
        const reports::AlignmentReportGenerator::ImageFileFormat statsImageFormat);

    /**
     * \brief Runs end-to-end alignment from the beginning
     */
    void run();


    enum State
    {
        Invalid = -2,
        Last = -1,
        Start = 0,          // constructor completed
        MatchFinderDone,    // MatchFinder done, foundMatchesMetadata_ is valid
        MatchSelectorDone,  // MatchSelector done, selectedMatchesMetadata_ is valid
        AlignmentReportsDone,
        BamDone,            // Bam file generated
        CasavaResetDone,    // Casava variant calling configured
        CasavaResumeDone,   // Casava variant calling completed
        Finish = CasavaResumeDone
    };



    AlignWorkflow::State getNextState() const;

    /**
     * \brief Performs single step of aligner state transition
     *
     * \return The new state
     */
    AlignWorkflow::State step();

    /**
     * \brief Changes the aligner state to the specified, provided the prerequisite data is available
     *
     * \return the new state
     */
    AlignWorkflow::State rewind(AlignWorkflow::State to);

private:
    template<class Archive> friend void serialize(Archive & ar, AlignWorkflow &, const unsigned int file_version);

    typedef alignment::BinMetadataList SelectedMatchesMetadata;

    const std::vector<std::string> &argv_;
    const std::string &casavaArgv_;
    const std::vector<flowcell::Layout> &flowcellLayoutList_;
    const bfs::path tempDirectory_;
    const bfs::path statsDirectory_;
    const bfs::path reportsDirectory_;
    const bfs::path projectsDirectory_;
    const bfs::path matchSelectorStatsXmlPath_;
    const unsigned coresMax_;
    const unsigned repeatThreshold_;
    const int mateDriftRange_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const bool allowVariableFastqLength_;
    const bool ignoreMissingBcls_;
    const bool ignoreMissingFilters_;
    const unsigned firstPassSeeds_;
    const unsigned long matchesPerBin_;
    const unsigned long availableMemory_;
    const unsigned mapqThreshold_;
    const bool pfOnly_;
    const unsigned baseQualityCutoff_;
    const bool keepUnaligned_;
    const bool putUnalignedInTheBack_;
    const bool clipSemialigned_;
    const unsigned gappedMismatchesMax_;
    const bool scatterRepeats_;
    const int gapMatchScore_;
    const int gapMismatchScore_;
    const int gapOpenScore_;
    const int gapExtendScore_;
    const alignment::TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore_;
    const unsigned inputLoadersMax_;
    const unsigned tempSaversMax_;
    const unsigned tempLoadersMax_;
    const unsigned outputSaversMax_;
    const build::GapRealignerMode realignGaps_;
    const int bamGzipLevel_;
    const double expectedBgzfCompressionRatio_;
    const bool keepDuplicates_;
    const std::string &binRegexString_;
    const common::ScoopedMallocBlock::Mode memoryControl_;
    const alignment::TemplateLengthStatistics userTemplateLengthStatistics_;
    const bfs::path demultiplexingStatsXmlPath_;
    const reports::AlignmentReportGenerator::ImageFileFormat statsImageFormat_;

    const reference::SortedReferenceXmlList sortedReferenceXmlList_;

    State state_;
    SelectedMatchesMetadata selectedMatchesMetadata_;
    build::BarcodeBamMapping barcodeBamMapping_;

    const alignWorkflow::FindMatchesTransition findMatchesTransition_;
    alignWorkflow::FoundMatchesMetadata foundMatchesMetadata_;

    static reference::SortedReferenceXmlList loadSortedReferenceXml(
        const reference::ReferenceMetadataList &referenceMetadataList);

    alignWorkflow::FoundMatchesMetadata findMatches() const;
    alignment::BinMetadataList selectMatches() const;
    void generateAlignmentReports() const;
    const build::BarcodeBamMapping generateBam(const SelectedMatchesMetadata &binPaths) const;
    void resetCasava(const build::BarcodeBamMapping & barcodeBamMapping) const;
    void resumeCasava(const build::BarcodeBamMapping & barcodeBamMapping) const;

};
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_HH
