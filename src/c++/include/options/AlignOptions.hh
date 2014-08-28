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
 ** \file AlignOptions.hh
 **
 ** Command line options for 'isaac-align'
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_OPTIONS_ALIGN_OPTIONS_HH
#define iSAAC_OPTIONS_ALIGN_OPTIONS_HH

#include <string>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#include "alignment/SeedMetadata.hh"
#include "build/GapRealigner.hh"
#include "common/Program.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "workflow/AlignWorkflow.hh"

namespace isaac
{
namespace options
{

class AlignOptions : public common::Options
{
public:
    AlignOptions();

    common::Options::Action parse(int argc, char *argv[]);

private:
    std::string usagePrefix() const {return "isaac-align";}
    void postProcess(boost::program_options::variables_map &vm);
    void verifyMandatoryPaths(boost::program_options::variables_map &vm);
    void parseParallelization();
    build::GapRealignerMode parseGapRealignment();
    void parseExecutionTargets();
    void parseMemoryControl();
    void parseGapScoring();
    workflow::AlignWorkflow::OptionalFeatures parseBamExcludeTags(std::string strBamExcludeTags);
    void parseDodgyAlignmentScore();
    void parseTemplateLength();
    void parseReferenceGenomes();
    std::vector<boost::filesystem::path> parseSampleSheetPaths() const;
    std::vector<std::pair<flowcell::Layout::Format, bool> > parseBaseCallsFormats();
    void parseStatsImageFormat();
    void parseQScoreBinValues();
    void parseBamExcludeTags();
    void processLegacyOptions(boost::program_options::variables_map &vm);

public:
    std::vector<std::string> argv;
    std::string description;
    std::vector<boost::filesystem::path> baseCallsDirectoryList;
    std::vector<std::string> baseCallsFormatStringList;
    std::vector<std::string> sampleSheetStringList;
    std::vector<std::string> barcodeMismatchesStringList;
    std::vector<std::string> tilesFilterList;
    std::vector<std::string> useBasesMaskList;
    std::vector<flowcell::Layout> flowcellLayoutList;
    flowcell::BarcodeMetadataList barcodeMetadataList;
    std::vector<boost::filesystem::path> sortedReferenceMetadataList;
    std::vector<std::string> referenceNameList;
    reference::ReferenceMetadataList referenceMetadataList;
    boost::filesystem::path tempDirectory;
    boost::filesystem::path outputDirectory;
    // the seed descriptor
    std::string seedDescriptor;
    unsigned seedLength;
    bool allowVariableFastqReadLength;
    bool allowVariableReadLength;
    unsigned laneNumberMax;
    bool cleanupIntermediary;
    bool ignoreMissingBcls;
    bool ignoreMissingFilters;
    // number of seeds to use on the first pass
    unsigned firstPassSeeds;
    // the list of seed metadata
    unsigned jobs;
    unsigned repeatThreshold;
    int mateDriftRange;
    unsigned neighborhoodSizeThreshold;
    std::string startFromString;
    workflow::AlignWorkflow::State startFrom;
    std::string stopAtString;
    workflow::AlignWorkflow::State stopAt;
    unsigned int verbosity;
    unsigned clustersAtATimeMax;
    bool ignoreNeighbors;
    bool ignoreRepeats;
    unsigned mapqThreshold;
    bool perTileTls;
    bool pfOnly;
    bool allowEmptyFlowcells_;
    unsigned baseQualityCutoff;
    std::string keepUnalignedString;
    bool keepUnaligned;
    bool preSortBins;
    bool putUnalignedInTheBack;
    bool realignGapsVigorously;
    bool realignDodgyFragments;
    unsigned realignedGapsPerFragment;
    bool clipSemialigned;
    bool clipOverlapping;
    bool scatterRepeats;
    unsigned gappedMismatchesMax;
    bool avoidSmithWaterman;
    std::string gapScoringString;
    int gapMatchScore;
    int gapMismatchScore;
    int gapOpenScore;
    int gapExtendScore;
    int minGapExtendScore;
    unsigned semialignedGapLimit;
    std::string dodgyAlignmentScoreString;
    alignment::TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore;
    std::string memoryControlString;
    common::ScoopedMallocBlock::Mode memoryControl;
    unsigned long memoryLimit;
    static const unsigned long memoryLimitUnlimited = 0;
    unsigned inputLoadersMax;
    unsigned tempSaversMax;
    unsigned tempLoadersMax;
    unsigned outputSaversMax;
    std::string realignGapsString;
    build::GapRealignerMode realignGaps;
    int bamGzipLevel;
    std::vector<std::string> bamHeaderTags;
    std::string bamPuFormat;
    double expectedBgzfCompressionRatio;
    bool singleLibrarySamples;
    bool keepDuplicates;
    bool markDuplicates;
    std::string binRegexString;
    std::vector<std::size_t> clusterIdList;
    alignment::TemplateLengthStatistics userTemplateLengthStatistics;
    std::string tlsString;
    std::vector<std::string> defaultAdapters;
    std::string statsImageFormatString;
    reports::AlignmentReportGenerator::ImageFileFormat statsImageFormat;
    bool bufferBins;
    bool qScoreBin;
    std::string qScoreBinValueString;
    boost::array<char, 256> fullBclQScoreTable;
    std::string bamExcludeTags;
    workflow::AlignWorkflow::OptionalFeatures optionalFeatures;
    bool pessimisticMapQ;
};

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_ALIGN_OPTIONS_HH
