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
    void parseGapRealignment();
    void parseExecutionTargets();
    void parseMemoryControl();
    void parseGapScoring();
    void parseKeepUnaligned();
    void parseDodgyAlignmentScore();
    void parseTemplateLength();
    void parseReferenceGenomes();
    std::vector<boost::filesystem::path> parseSampleSheetPaths() const;
    std::vector<flowcell::Layout::Format> parseBaseCallsFormats();
    const std::vector<const char*> filterCasavaOptions(char * const *begin, char * const *end);
    void parseStatsImageFormat();

public:
    std::vector<std::string> argv;
    std::string casavaArgv;
    std::vector<boost::filesystem::path> baseCallsDirectoryList;
    std::vector<std::string> baseCallsFormatStringList;
    std::vector<std::string> sampleSheetStringList;
    std::vector<std::string> barcodeMismatchesStringList;
    std::vector<std::string> tilesFilterList;
    std::vector<std::string> useBasesMaskList;
    std::vector<flowcell::Layout> flowcellLayoutList;
    flowcell::BarcodeMetadataList barcodeMetadataList;
    std::vector<boost::filesystem::path> sortedReferenceXmlList;
    std::vector<std::string> referenceNameList;
    reference::ReferenceMetadataList referenceMetadataList;
    boost::filesystem::path tempDirectory;
    boost::filesystem::path outputDirectory;
    // the seed descriptor
    std::string seedDescriptor;
    bool allowVariableFastqLength;
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
    bool ignoreNeighbors;
    bool ignoreRepeats;
    unsigned mapqThreshold;
    bool pfOnly;
    unsigned baseQualityCutoff;
    std::string keepUnalignedString;
    bool keepUnaligned;
    bool putUnalignedInTheBack;
    bool clipSemialigned;
    unsigned gappedMismatchesMax;
    bool scatterRepeats;
    std::string gapScoringString;
    int gapMatchScore;
    int gapMismatchScore;
    int gapOpenScore;
    int gapExtendScore;
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
    double expectedBgzfCompressionRatio;
    bool keepDuplicates;
    std::string binRegexString;
    std::vector<size_t> clusterIdList;
    alignment::TemplateLengthStatistics userTemplateLengthStatistics;
    std::string tlsString;
    std::vector<std::string> defaultAdapters;
    std::string statsImageFormatString;
    reports::AlignmentReportGenerator::ImageFileFormat statsImageFormat;

};

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_ALIGN_OPTIONS_HH
