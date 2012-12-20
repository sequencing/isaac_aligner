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
 ** \file AlignOptions.cpp
 **
 ** Command line options for 'isaac'
 **
 ** \author Come Raczy
 **/

#include <algorithm>
#include <string>
#include <vector>
#include <ostream>
#include <fstream>

#include <boost/assign.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/algorithm/string/regex.hpp>

#include "basecalls/ConfigXml.hh"
#include "common/Exceptions.hh"
#include "demultiplexing/SampleSheetCsv.hh"
#include "oligo/Mask.hh"
#include "options/AlignOptions.hh"

#include "alignOptions/BclFlowcell.hh"
#include "alignOptions/FastqFlowcell.hh"
#include "alignOptions/DefaultAdaptersOption.hh"

namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;
using common::InvalidOptionException;
using boost::format;

/**
 * \brief Throws an exception in case of troubles,
 *
 * \return returns ulimit -v in bytes or 0 if unlimited
 */
unsigned long getUlimitV()
{
    unsigned long ret = 0;
    if (!isaac::common::ulimitV(&ret))
    {
        BOOST_THROW_EXCEPTION(isaac::common::ResourceException(
            errno, (boost::format("Failed to get the memory consumption limit: %s") % strerror(errno)).str()));
    }
    return ret == (0UL -1) ? 0 : ret;
}

AlignOptions::AlignOptions()
    : baseCallsDirectoryList()
    , baseCallsFormatStringList(1, "bcl")
    , barcodeMismatchesStringList(1, "0")
    , referenceNameList(1, "default")
    , tempDirectory("./Temp")
    , outputDirectory("./Aligned")
    , seedDescriptor("16:0:32:64")
    , allowVariableFastqLength(false)
    , ignoreMissingBcls(false)
    , ignoreMissingFilters(false)
    , firstPassSeeds(1)
    , jobs(boost::thread::hardware_concurrency())
    , repeatThreshold(10)
    , mateDriftRange(-1)
    , neighborhoodSizeThreshold(0) //(10000) - disabled by default as so far the neighbor matcher only increased the probability of misplaced reads
    , startFromString("Start")
    , startFrom(workflow::AlignWorkflow::Start)
    , stopAtString("Finish")
    , stopAt(workflow::AlignWorkflow::Finish)
    , verbosity(2)
    , ignoreNeighbors(false)
    , ignoreRepeats(false)
    , mapqThreshold(0)
    , pfOnly(true)
    , baseQualityCutoff(20)
    , keepUnalignedString("discard")
    , keepUnaligned(false)
    , putUnalignedInTheBack(false)
    , clipSemialigned(true)
    , gappedMismatchesMax(5)
    , scatterRepeats(false)
    , gapScoringString("eland")
    , gapMatchScore(2)
    , gapMismatchScore(-1)
    , gapOpenScore(-15)
    , gapExtendScore(-3)
    , dodgyAlignmentScoreString("Zero")
    , dodgyAlignmentScore(alignment::TemplateBuilder::Zero)
#ifdef ISAAC_THREAD_CERR_DEV_TRACE_ENABLED
    , memoryControlString("off")
#else //ISAAC_THREAD_CERR_DEV_TRACE_ENABLED
    , memoryControlString("warning")
#endif //ISAAC_THREAD_CERR_DEV_TRACE_ENABLED
    , memoryControl(common::ScoopedMallocBlock::Invalid)
    , memoryLimit(getUlimitV() / 1024 / 1024 / 1024)
    , inputLoadersMax(64) // bcl files are small, there are lots of them and at the moment they are expected to sit on a highly-parallelizable high-latency network storage
    , tempSaversMax(64)   // currently most runs using isilon as Temp storage. In this case fragmentation is not an issue
    , tempLoadersMax(8)
    , outputSaversMax(8)
    , realignGapsString("yes")
    , bamGzipLevel(boost::iostreams::gzip::best_speed)
    , expectedBgzfCompressionRatio(1)
    , keepDuplicates(false)
    , binRegexString()
    , userTemplateLengthStatistics(-1)
    , statsImageFormatString("gif")
    , statsImageFormat(reports::AlignmentReportGenerator::gif)
{
    namedOptions_.add_options()
        ("base-calls-directory,b"   , bpo::value<std::vector<bfs::path> >(&baseCallsDirectoryList)->multitoken(),
                "full path to the base calls directory. Multiple entries allowed.")
        ("base-calls-format"        , bpo::value<std::vector<std::string> >(&baseCallsFormatStringList)->multitoken(),
                "Multiple entries allowed. Each entry is applied to the corresponding base-calls-directory. "
                "Last entry is applied to all --base-calls-directory that don't have --base-calls-format specified."
                "\n  - bcl(default)    : common bcl files, no compression."
                "\n  - bcl-gz          : bcl files are individually compressed and named s_X_YYYY.bcl.gz"
                "\n  - fastq           : One fastq per lane/read named lane<X>_read<Y>.fastq and located directly in "
                                        "the specified base-calls-directory"
                "\n  - fastq-gz        : One compressed fastq per lane/read named lane<X>_read<Y>.fastq.gz and "
                                        "located directly in the specified base-calls-directory")
        ("default-adapters"
                                    , bpo::value<std::vector<std::string> >(&defaultAdapters)->multitoken(),
                "Multiple entries allowed. Each entry is associated with the corresponding base-calls-directory. "
                "Flowcells that don't have default-adapters provided, don't get adapters clipped in the data. "
                "\nEach entry is a comma-separated list of adapter sequences written in the direction of the reference. "
                "Wildcard (* character) is allowed only on one side of the sequence. "
                "Entries with * apply only to the alignments on the matching strand. Entries without * apply to all "
                "strand alignments and are matched in the order of appearance in the list."
                "\nExamples:"
                "\n  ACGT*,*TGCA       : Will clip ACGT and all subsequent bases in the forward-strand alignments "
                                       "and mirror the behavior for the reverse-strand alignments."
                "\n  ACGT,TGCA         : Will find the following sequences in the reads: ACGT, TGCA, ACGTTGCA "
                                        " (but not TGCAACGT!) regardless "
                                        "of the alignment strand. Then will attempt to clip off the side of the read "
                                        "that is shorter. If both sides are roughly equal length, will clip off the "
                                        "side that has less matches."
                "\n  Nextera           : Nextera standard, Same as CTGTCTCTTATACACATCT*,*AGATGTGTATAAGAGACAG"
                "\n  NexteraMp         : Nextera mate-pair. Same as CTGTCTCTTATACACATCT,AGATGTGTATAAGAGACAG")
        ("sample-sheet,s"   , bpo::value<std::vector<std::string> >(&sampleSheetStringList)->multitoken(),
                "Multiple entries allowed. Each entry is applied to the corresponding base-calls-directory."
                "\n  - none            : process flowcell as if there is no sample sheet"
                "\n  - default         : use <base-calls-directory>/SampleSheet.csv if it exists. This is the default behavior."
                "\n  - <file path>     : use <file path> as sample sheet for the flowcell.")
        ("barcode-mismatches"   , bpo::value<std::vector<std::string> >(&barcodeMismatchesStringList)->multitoken()->default_value(barcodeMismatchesStringList, barcodeMismatchesStringList.at(0)),
                "Multiple entries allowed. Each entry is applied to the corresponding base-calls-directory. "
                "Last entry applies to all the bases-calls-directory that do not have barcode-mismatches specified. "
                "Last component mismatch value applies to all subsequent barcode components should there be more than one. "
                "Examples: "
                "\n  - 1:0             : allow one mismatch for the first barcode component and no mismatches for the subsequent components."
                "\n  - 1               : allow one mismatch for every barcode component."
                "\n  - 0               : no mismatches allowed in any barcode component. This is the default.")
        ("realign-gaps"         , bpo::value<std::string>(&realignGapsString)->default_value(realignGapsString),
                "For reads overlapping the gaps occurring on other reads, check if applying "
                "those gaps reduces mismatch count. Significantly reduces number of false SNPs "
                "reported around short indels."
                "\n  - no              : no gap realignment"
                "\n  - yes             : discard all existing gaps in the read and apply the combination of gaps found in the region that "
                "yields maximum reduction of mismatches")
        ("bam-gzip-level"           , bpo::value<int>(&bamGzipLevel)->default_value(bamGzipLevel),
                "Gzip level to use for BAM")
        ("expected-bgzf-ratio"           , bpo::value<double>(&expectedBgzfCompressionRatio)->default_value(expectedBgzfCompressionRatio),
                "compressed = ratio * uncompressed. To avoid memory overallocation during the bam generation, iSAAC has to assume certain compression ratio. "
                "If iSAAC estimates less memory than is actually required, it will fail at runtime. You can check how far "
                "you are from the dangerous zone by looking at the resident/swap memory numbers for your process "
                "during the bam generation. If you see too much showing as 'swap', it is safe to reduce the --expected-bgzf-ratio.")
        ("tiles"                    , bpo::value<std::vector<std::string> >(&tilesFilterList),
                "Comma-separated list of regular expressions to select only a subset of the tiles available in the flow-cell."
                "\n- to select all the tiles ending with '5' in all lanes: --tiles [0-9][0-9][0-9]5"
                "\n- to select tile 2 in lane 1 and all the tiles in the other lanes: --tiles s_1_0002,s_[2-8]"
                "\nMultiple entries allowed, each applies to the corresponding base-calls-directory.")
        ("use-bases-mask"           , bpo::value<std::vector<std::string> >(&useBasesMaskList),
                "Conversion mask characters:"
                "\n  - Y or y          : use"
                "\n  - N or n          : discard"
                "\n  - I or i          : use for indexing"
                "\n"
                "\nIf not given, the mask will be guessed from the B<config.xml> file in the base-calls directory."
                "\n"
                "\nFor instance, in a 2x76 indexed paired end run, the mask I<Y76,I6n,y75n> means:"
                "\n  use all 76 bases from the first end, discard the last base of the indexing read, and "
                "use only the first 75 bases of the second end.")
        ("seeds"                  , bpo::value<std::string>(&seedDescriptor)->default_value(seedDescriptor),
                "Seed descriptors for each read, given as a comma-separated list-of-seeds for each read. "
                "A list-of-seeds is a colon-separated list of offsets from the beginning of the read. "
                "\nExamples:"
                "\n  - 0:32,0:32:64    : two seeds on the first read (at offsets 0 and 32) "
                "and three seeds on the second read (at offsets 0, 32, and 64) and on subsequent reads."
                "\n  - 0:32:64         : three seeds on all the reads (at offsets 0, 32 and 64)"
                "\nNote that the last list-of-seeds is repeated to all subsequent reads "
                "if there are more reads than there are colon-separated lists-of-seeds."
            )
        ("first-pass-seeds,f"       , bpo::value<unsigned>(&firstPassSeeds)->default_value(firstPassSeeds),
                "the number of seeds to use in the first pass of the match finder")
        ("reference-genome,r"       , bpo::value<std::vector<bfs::path> >(&sortedReferenceXmlList),
                "Full path to the reference genome XML descriptor. Multiple entries allowed."
                "Each entry applies to the corresponding --reference-name. The last --reference-genome entry "
                "may not have a corresponding --reference-name. In this case the default name 'default' is assumed."
            )
        ("reference-name,n"       , bpo::value<std::vector<std::string> >(&referenceNameList),
                "Unique symbolic name of the reference. Multiple entries allowed. Each entry is associated with "
                "the corresponding --reference-genome and will be matched against the 'reference' column "
                "in the sample sheet. "
                "\nSpecial names:"
                "\n  - unknown         : default reference to use with data that did not match any barcode."
                "\n  - default         : reference to use for the data with no matching value in sample sheet 'reference' column."
            )
        ("temp-directory,t"         , bpo::value<bfs::path>(&tempDirectory)->default_value(tempDirectory),
                "Directory where the temporary files will be stored (matches, unsorted alignments, etc.)")
        ("output-directory,o"       , bpo::value<bfs::path>(&outputDirectory)->default_value(outputDirectory),
                "Directory where the final alignment data be stored")
        ("jobs,j"                   , bpo::value<unsigned int>(&jobs)->default_value(jobs),
                "Maximum number of compute threads to run in parallel")
        ("input-parallel-load"            , bpo::value<unsigned>(&inputLoadersMax)->default_value(inputLoadersMax),
                "Maximum number of parallel file read operations for --base-calls-directory")
        ("temp-parallel-load"            , bpo::value<unsigned>(&tempLoadersMax)->default_value(tempLoadersMax),
                "Maximum number of parallel file read operations for --temp-directory")
        ("temp-parallel-save"            , bpo::value<unsigned>(&tempSaversMax)->default_value(tempSaversMax),
                "Maximum number of parallel file write operations for --temp-directory")
        ("output-parallel-save"            , bpo::value<unsigned>(&outputSaversMax)->default_value(outputSaversMax),
                "Maximum number of parallel file write operations for --output-directory")
        ("repeat-threshold"         , bpo::value<unsigned>(&repeatThreshold)->default_value(repeatThreshold),
                "Threshold used to decide if matches must be discarded as too abundant (when the number of repeats is greater or equal to the threshold)")
        ("shadow-scan-range"         , bpo::value(&mateDriftRange)->default_value(mateDriftRange),
                "-1     - scan for possible mate alignments between template min and max\n"
                ">=0    - scan for possible mate alignments in range of template median += shadow-scan-range")
        ("neighborhood-size-threshold", bpo::value<unsigned>(&neighborhoodSizeThreshold)->default_value(neighborhoodSizeThreshold),
                "Threshold used to decide if the number of reference 32-mers sharing the same prefix (16 bases) "
                "is small enough to justify the neighborhood search. Use large enough value e.g. 10000 to enable "
                "alignment to positions where seeds don't match exactly.")
        ("verbosity"                , bpo::value<unsigned int>(&verbosity)->default_value(verbosity),
                "Verbosity: FATAL(0), ERRORS(1), WARNINGS(2), INFO(3), DEBUG(4) (not supported yet)")
        ("start-from"               , bpo::value<std::string>(&startFromString)->default_value(startFromString),
                "Start processing at the specified stage:"
                "\n  - Start            : don't resume, start from beginning"
                "\n  - MatchFinder      : same as Start"
                "\n  - MatchSelector    : skip match identification, continue with template selection"
                "\n  - AlignmentReports : regenerate alignment reports and bam"
                "\n  - Bam              : resume at bam generation"
#ifdef HAVE_CASAVA
                "\n  - CasavaReset      : resume at CASAVA configuration (will reset CASAVA analyses)"
                "\n  - CasavaResume     : resume at CASAVA execution (will resume interrupted CASAVA execution)"
#endif
                "\n  - Finish           : Same as CasavaResume."
                "\n  - Last             : resume from the last successful step"
                "\nNote that although iSAAC attempts to perform some basic validation, the only safe option is 'Start' "
                "The primary purpose of the feature is to reduce the time required to diagnose the issues rather than "
                "be used on a regular basis."
        )
        ("stop-at"               , bpo::value<std::string>(&stopAtString)->default_value(stopAtString),
                "Stop processing after the specified stage is complete:"
                "\n  - Start            : perform the first stage only"
                "\n  - MatchFinder      : same as Start"
                "\n  - MatchSelector    : don't generate alignment reports and bam"
                "\n  - AlignmentReports : don't perform bam generation"
                "\n  - Bam              : finish when bam is done"
#ifdef HAVE_CASAVA
                "\n  - CasavaReset      : finish after CASAVA configuration (will reset CASAVA analyses)"
                "\n  - CasavaResume     : finish after CASAVA execution completes"
#endif
                "\n  - Finish           : stop at the end."
                "\n  - Last             : perform up to the last successful step only"
                "\nNote that although iSAAC attempts to perform some basic validation, the only safe option is 'Finish' "
                "The primary purpose of the feature is to reduce the time required to diagnose the issues rather than "
                "be used on a regular basis."
        )
        ("ignore-neighbors"         , bpo::value<bool>(&ignoreNeighbors)->default_value(ignoreNeighbors),
                "When not set, MatchFinder will ignore perfect seed matches during single-seed pass, "
                "if the reference k-mer is known to have neighbors.")
        ("ignore-repeats"         , bpo::value<bool>(&ignoreRepeats)->default_value(ignoreRepeats),
                "Normally exact repeat matches prevent inexact seed matching. If this flag is set, inexact "
                "matches will be considered even for the seeds that match to repeats.")
        ("mapq-threshold"           , bpo::value<unsigned>(&mapqThreshold)->default_value(mapqThreshold),
                "Threshold used to filter the templates based on their mapping quality: the BAM file will only "
                "contain the templates with a mapping quality greater than or equal to the threshold. Templates "
                "(or fragments) with a mapping quality of 4 or more are guaranteed to be uniquely aligned. "
                "Those with a mapping quality of 3 or less are either mapping to repeat regions or have a large number of errors.")
        ("pf-only"                  , bpo::value<bool>(&pfOnly)->default_value(pfOnly),
                "When set, only the fragments passing filter (PF) are generated in the BAM file")
        ("scatter-repeats"                  , bpo::value<bool>(&scatterRepeats)->default_value(scatterRepeats),
                "When set, extra care will be taken to scatter pairs aligning to repeats across the repeat locations ")
        ("base-quality-cutoff"     , bpo::value<unsigned>(&baseQualityCutoff)->default_value(baseQualityCutoff),
                "3' end quality trimming cutoff. 0 turns off the trimming.")
        ("variable-fastq-read-length"  , bpo::value<bool>(&allowVariableFastqLength)->default_value(allowVariableFastqLength),
                "Unless set, iSAAC will fail if the length of the sequence changes between the records of the fastq file.")
        ("ignore-missing-bcls"      , bpo::value<bool>(&ignoreMissingBcls)->default_value(ignoreMissingBcls),
                "When set, missing bcl files are treated as all clusters having N bases for the "
                "corresponding tile cycle. Otherwise, encountering a missing bcl file causes the analysis to fail.")
        ("ignore-missing-filters"      , bpo::value<bool>(&ignoreMissingFilters)->default_value(ignoreMissingFilters),
                "When set, missing filter files are treated as if all clusters pass filter for the "
                "corresponding tile. Otherwise, encountering a missing filter file causes the analysis to fail.")
        ("keep-unaligned"           , bpo::value<std::string>(&keepUnalignedString)->default_value(keepUnalignedString),
                "Available options:"
                "\n - discard          : discard clusters where both reads are not aligned"
                "\n - front            : keep unaligned clusters in the front of the BAM file"
                "\n - back             : keep unaligned clusters in the back of the BAM file")
        ("clip-semialigned"         , bpo::value<bool>(&clipSemialigned)->default_value(clipSemialigned),
                "When set, semialigned reads have their mismatching part soft-clipped")
        ("gapped-mismatches"   , bpo::value<unsigned>(&gappedMismatchesMax)->default_value(gappedMismatchesMax),
                "Maximum number of mismatches allowed to accept a gapped alignment.")
        ("gap-scoring"   , bpo::value<std::string>(&gapScoringString)->default_value(gapScoringString),
                "Gapped alignment algorithm parameters:"
                "\n - eland            : equivalent of 2:-1:-15:-3"
                "\n - bwa              : equivalent of 0:-3:-11:-4"
                "\n - m:mm:go:ge       : colon-delimited string of values where:"
                "\n     m              : match score"
                "\n     mm             : mismatch score"
                "\n     go             : gap open score"
                "\n     ge             : gap extend score")
        ("dodgy-alignment-score"   , bpo::value<std::string>(&dodgyAlignmentScoreString)->default_value(dodgyAlignmentScoreString),
                "Controls the behavior for templates where alignment score is impossible to assign:"
                "\n - Unaligned        : marks template fragments as unaligned"
                "\n - Zero             : sets the fragment and pair alignment scores to 0"
                "\n - Unknown          : assigns value 255 for bam MAPQ. Ensures SM and AS are not specified in the bam"
        )
        ("keep-duplicates" , bpo::value<bool>(&keepDuplicates)->default_value(keepDuplicates),
                "Keep duplicate pairs in the bam file (with 0x400 flag set in all but the best one)")
        ("bin-regex" , bpo::value<std::string>(&binRegexString)->default_value(binRegexString),
                "Comma-separated list of regular expressions. If not empty, bam files will be filtered to contain "
                "only the bins that match.")
        ("memory-control"           , bpo::value<std::string>(&memoryControlString)->default_value(memoryControlString),
                "Define the behavior in case unexpected memory allocations are detected: "
                "\n  - warning         : Log WARNING about the allocation."
                "\n  - off             : Don't monitor dynamic memory usage."
                "\n  - strict          : Fail memory allocation. Intended for development use."
        )
        ("memory-limit,m"           , bpo::value<unsigned long>(&memoryLimit)->default_value(memoryLimit),
                "Limits major memory consumption operations to a set number of gigabytes. "
                "0 means no limit, however 0 is not allowed as in such case iSAAC will most likely consume "
                "all the memory on the system and cause it to crash. Default value is taken from ulimit -v.")
        ("cluster,c"                , bpo::value<std::vector<size_t> >(&clusterIdList)->multitoken(),
                "Restrict the alignment to the specified cluster Id (multiple entries allowed)")
        ("tls"                      , bpo::value<std::string>(&tlsString),
                "Template-length statistics in the format 'min:median:max:lowStdDev:highStdDev:M0:M1', "
                "where M0 and M1 are the numeric value of the models (0=FFp, 1=FRp, 2=RFp, 3=RRp, 4=FFm, 5=FRm, 6=RFm, 7=RRm)")
#ifdef HAVE_CASAVA
        ("casava*", "Any option beginning with --casava will have the --casava prefix removed and passed along "
                    "with its arguments to the configureBuild.pl CASAVA script at the CasavaReset stage. "
                    "For example --casava--variantsConsensusVCF will add --variantsConsensusVCF to the configureBuild.pl "
                    "command line")
#endif
        ("stats-image-format", bpo::value<std::string>(&statsImageFormatString)->default_value(statsImageFormatString),
                "Format to use for images during stats generation"
                "\n - gif        : produce .gif type plots"
                "\n - none       : no stat generation"
        )
        ;



}

/**
 * \brief Extracts the number of allowed mismatches for each barcode component from command line arg.
 *
 */
static std::vector<unsigned> parseBarcodeMismatches(
    const std::string &barcodeMismatches)
{
    std::vector<std::string> barcodeMismatchesList;
    boost::algorithm::split(barcodeMismatchesList, barcodeMismatches,  boost::algorithm::is_any_of(":"));
    std::vector<unsigned> ret;
    std::transform(barcodeMismatchesList.begin(), barcodeMismatchesList.end(),
                   std::back_inserter(ret),
                   &boost::lexical_cast<unsigned, std::string>);

    return ret;
}

void AlignOptions::verifyMandatoryPaths(bpo::variables_map &vm)
{
    const std::vector<std::string> requiredOptions = boost::assign::list_of("reference-genome")("base-calls-directory");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const format message = format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
    typedef std::pair<bfs::path *, std::string> PathOption;
    BOOST_FOREACH(bfs::path &baseCallsDirectory, baseCallsDirectoryList)
    {
        if (baseCallsDirectory.empty())
        {
            const format message = format("\n   *** The 'base-calls-directory' can't be empty (use '.' for current directory) ***\n") % baseCallsDirectory;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        if(!exists(baseCallsDirectory))
        {
            const format message = format("\n   *** The 'base-calls-directory' does not exist: %s ***\n") % baseCallsDirectory;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        baseCallsDirectory = boost::filesystem::absolute(baseCallsDirectory);
    }

    BOOST_FOREACH(bfs::path &sortedReferenceXml, sortedReferenceXmlList)
    {
        if (sortedReferenceXml.empty())
        {
            const format message = format("\n   *** The 'reference-genome' can't be empty ***\n") % sortedReferenceXml;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        if(!exists(sortedReferenceXml))
        {
            const format message = format("\n   *** The 'reference-genome' does not exist: %s ***\n") % sortedReferenceXml;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        sortedReferenceXml = boost::filesystem::absolute(sortedReferenceXml);
        if(!exists(sortedReferenceXml))
        {
            const format message = format("\n   *** The 'reference-genome' does not exist: %s ***\n") % sortedReferenceXml;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }

    const std::vector<PathOption> pathOptions = boost::assign::list_of
        (PathOption(&tempDirectory, "temp-directory"))
        (PathOption(&outputDirectory, "output-directory"));
    BOOST_FOREACH(const PathOption &pathOption, pathOptions)
    {
        if(pathOption.first->empty())
        {
            const format message = format("\n   *** The '%s' can't be empty (use '.' for current directory) ***\n") % pathOption.second;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        *pathOption.first = boost::filesystem::absolute(*pathOption.first);
    }
}

void AlignOptions::parseParallelization()
{
    if(0 >= jobs)
    {
        const format message = format("\n   *** The 'jobs' option must be strictly positive (%d) ***\n") % jobs;
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
}

void AlignOptions::parseGapRealignment()
{
    if(realignGapsString == "no")
    {
        realignGaps = build::GAP_REALIGNER_OFF;
    }
    else if(realignGapsString == "yes")
    {
        realignGaps = build::GAP_REALIGNER_ON;
    }
    else
    {
        const format message = format("\n   *** The 'realign-gaps' value is invalid %s ***\n") % realignGapsString;
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
}

void AlignOptions::parseExecutionTargets()
{
    const static std::vector<std::string> allowedStageStrings =
        boost::assign::list_of("Start")("MatchFinder")("MatchSelector")
            ("AlignmentReports")("Bam")
            ("CasavaReset")("CasavaResume")
            ("Finish")("Last");

#ifndef HAVE_CASAVA
    if (startFromString == "CasavaReset" || startFromString == "CasavaResume")
    {
        const boost::format message = boost::format("\n   *** Invalid value given '%s' for --start-from. CASAVA integration is not enabled for this installation. ***\n") %
                startFromString;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
    if (stopAtString == "CasavaReset" || stopAtString == "CasavaResume")
    {
        const boost::format message = boost::format("\n   *** Invalid value given '%s' for --stop-at. CASAVA integration is not enabled for this installation. ***\n") %
            stopAtString;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
#endif

    std::vector<std::string>::const_iterator startFromIt =
        std::find(allowedStageStrings.begin(), allowedStageStrings.end(), startFromString);
    if (allowedStageStrings.end() == startFromIt)
    {
        const boost::format message = boost::format("\n   *** Invalid value given '%s' for --start-from ***\n") %
                startFromString;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    const std::vector<std::string>::const_iterator::difference_type startFromPos =
        startFromIt - allowedStageStrings.begin();
    startFrom =
        0 == startFromPos ? workflow::AlignWorkflow::Start :
        1 == startFromPos ? workflow::AlignWorkflow::Start :
        2 == startFromPos ? workflow::AlignWorkflow::MatchFinderDone :
        3 == startFromPos ? workflow::AlignWorkflow::MatchSelectorDone :
        4 == startFromPos ? workflow::AlignWorkflow::AlignmentReportsDone :
        5 == startFromPos ? workflow::AlignWorkflow::BamDone :
        6 == startFromPos ? workflow::AlignWorkflow::CasavaResetDone :
        7 == startFromPos ? workflow::AlignWorkflow::CasavaResetDone :
                            workflow::AlignWorkflow::Last;

    std::vector<std::string>::const_iterator stopAtIt =
        std::find(allowedStageStrings.begin(), allowedStageStrings.end(), stopAtString);
    if (allowedStageStrings.end() == stopAtIt)
    {
        const boost::format message = boost::format("\n   *** Invalid value given '%s' for --stop-at ***\n") %
            stopAtString;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    const std::vector<std::string>::const_iterator::difference_type stopAtPos =
        stopAtIt - allowedStageStrings.begin();
    stopAt =
        0 == stopAtPos ? workflow::AlignWorkflow::MatchFinderDone :
        1 == stopAtPos ? workflow::AlignWorkflow::MatchFinderDone :
        2 == stopAtPos ? workflow::AlignWorkflow::MatchSelectorDone :
        3 == stopAtPos ? workflow::AlignWorkflow::AlignmentReportsDone :
        4 == stopAtPos ? workflow::AlignWorkflow::BamDone :
        5 == stopAtPos ? workflow::AlignWorkflow::CasavaResetDone :
        6 == stopAtPos ? workflow::AlignWorkflow::Finish :
        7 == stopAtPos ? workflow::AlignWorkflow::Finish :
                         workflow::AlignWorkflow::Last;
}

void AlignOptions::parseMemoryControl()
{
    const std::vector<std::string> allowedMemoryControlStrings =
        boost::assign::list_of("warning")("strict")("off");
    std::vector<std::string>::const_iterator memoryControlIt =
        std::find(allowedMemoryControlStrings.begin(), allowedMemoryControlStrings.end(), memoryControlString);
    if (allowedMemoryControlStrings.end() == memoryControlIt)
    {
        const boost::format message = boost::format("\n   *** Invalid value given '%s' for --memory-control ***\n") %
            memoryControlString;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
    const std::vector<std::string>::const_iterator::difference_type memoryControlPos =
        memoryControlIt - allowedMemoryControlStrings.begin();
    memoryControl =
        0 == memoryControlPos ? common::ScoopedMallocBlock::Warning :
        1 == memoryControlPos ? common::ScoopedMallocBlock::Strict : common::ScoopedMallocBlock::Off;

    if (!memoryLimit)
    {
        const format message = format("\n   *** The 'memory-limit' option must be a positive value in gigabytes ***\n");
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
}

void AlignOptions::parseGapScoring()
{
    if ("bwa" == gapScoringString)
    {
        gapMatchScore = 0;
        gapMismatchScore = -3;
        gapOpenScore = -11;
        gapExtendScore = -4;
    }
    else if ("eland" == gapScoringString)
    {
        gapMatchScore = 2;
        gapMismatchScore = -1;
        gapOpenScore = -15;
        gapExtendScore = -3;
    }
    else
    {
        std::vector<std::string> gapScores;
        using boost::algorithm::split;
        using boost::algorithm::is_any_of;

        boost::split(gapScores, gapScoringString,  is_any_of(":"));

        if (4 != gapScores.size())
        {
            const format message = format("\n   *** The 'gap-scoring' string must contain four components delimited by ':' ***\n");
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }

        gapMatchScore = boost::lexical_cast<int>(gapScores.at(0));
        gapMismatchScore = boost::lexical_cast<int>(gapScores.at(1));
        gapOpenScore = boost::lexical_cast<int>(gapScores.at(2));
        gapExtendScore = boost::lexical_cast<int>(gapScores.at(3));
    }
}

void AlignOptions::parseKeepUnaligned()
{
    if ("discard" == keepUnalignedString)
    {
        keepUnaligned = false;
    }
    else if ("front" == keepUnalignedString)
    {
        keepUnaligned = true;
        putUnalignedInTheBack = false;
    }
    else if ("back" == keepUnalignedString)
    {
        keepUnaligned = true;
        putUnalignedInTheBack = true;
    }
    else
    {
        const format message = format("\n   *** The 'keep-unaligned' string must must be 'discard', 'front' or 'back'***\n");
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
}

void AlignOptions::parseDodgyAlignmentScore()
{
    if ("Unknown" == dodgyAlignmentScoreString)
    {
        dodgyAlignmentScore = alignment::TemplateBuilder::Unknown;
    }
    else if ("Zero" == dodgyAlignmentScoreString)
    {
        dodgyAlignmentScore = alignment::TemplateBuilder::Zero;
    }
    else if ("Unaligned" == dodgyAlignmentScoreString)
    {
        dodgyAlignmentScore = alignment::TemplateBuilder::Unaligned;
    }
    else
    {
        const format message = format("\n   *** The 'dodgy-alignment-score' option must be either Unknown, Unaligned or Zero (%s given) ***\n") % dodgyAlignmentScoreString;
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
}
void AlignOptions::parseTemplateLength()
{
    if (!tlsString.empty())
    {
        std::vector<std::string> tlsParameters;
        boost::split_regex(tlsParameters, tlsString, boost::regex(":"));
        if (7 != tlsParameters.size())
        {
            const format message = format("\n   *** The 'tls' option must be 'min:median:max:low:high:M0:M1' (%s) ***\n") % tlsString;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        std::vector<unsigned> v;
        v.reserve(tlsParameters.size());
        BOOST_FOREACH(const std::string &p, tlsParameters)
        {
            try
            {
                v.push_back(boost::lexical_cast<unsigned>(p));
            }
            catch(...)
            {
                const format message = format("\n   *** All components of the 'tls' must be unsigned integers: found '%s' in '%s' ***\n") % p % tlsString;
                BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
            }
        }
        using alignment::TemplateLengthStatistics;
        if (std::max(v[5], v[6]) >= TemplateLengthStatistics::InvalidAlignmentModel)
        {
            const format message = format("\n   *** Both models in the 'tls' must be less than %d (%s) ***\n") %
                TemplateLengthStatistics::InvalidAlignmentModel % tlsString;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        typedef TemplateLengthStatistics::AlignmentModel AlignmentModel;
        const AlignmentModel M0 = static_cast<AlignmentModel>(v[5]);
        const AlignmentModel M1 = static_cast<AlignmentModel>(v[6]);
        userTemplateLengthStatistics = TemplateLengthStatistics(v[0], v[2], v[1], v[3], v[4], M0, M1, true);
    }
}

void AlignOptions::parseReferenceGenomes()
{
    if (referenceNameList.size() > sortedReferenceXmlList.size())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Too many --reference-name options specified. There must be at most one per --reference-genome. ***\n"));
    }
    std::vector<std::string> dupeChecking(referenceNameList);
    std::sort(dupeChecking.begin(), dupeChecking.end());
    std::vector<std::string>::const_iterator firstDupe =
        std::unique(dupeChecking.begin(), dupeChecking.end());

    if (dupeChecking.end() != firstDupe)
    {
        std::string msg =
            (boost::format("\n   *** --reference-name arguments must be unique. '%s' is supplied more than once. ***\n")
            % *firstDupe).str();
        BOOST_THROW_EXCEPTION(InvalidOptionException(msg));
    }

    if (referenceNameList.size() < sortedReferenceXmlList.size())
    {
        if (referenceNameList.end() == std::find(referenceNameList.begin(), referenceNameList.end(), "default"))
        {
            referenceNameList.push_back("default");
        }
    }

    if (referenceNameList.size() < sortedReferenceXmlList.size())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Not enough --reference-name options specified. There must be one for each --reference-genome. ***\n"));
    }

    BOOST_FOREACH(const boost::filesystem::path &sortedReferenceXml, sortedReferenceXmlList)
    {
        const unsigned referenceIndex = &sortedReferenceXml - &sortedReferenceXmlList.front();
        referenceMetadataList.push_back(
            reference::ReferenceMetadata(referenceNameList.at(referenceIndex), sortedReferenceXml,
                                         referenceMetadataList.size()));
    }
}

/**
 * \brief extracts all command line parameters that begin with --casava and their argument
 *
 * \return returns vector of everything but the options beginning --casava and their arguments
 */
const std::vector<const char*> AlignOptions::filterCasavaOptions(char * const*begin, char *const*end)
{
    std::vector<const char*> ret;

    bool extractCasavaArgs = false;
    BOOST_FOREACH(const char *arg, std::make_pair(begin, end))
    {
        if (!strncmp("--casava*", arg, 9))
        {
            BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Invalid use of --casava* option. Please read the option help again. ***\n"));
        }

        if (!strncmp("--casava", arg, 8))
        {
            arg += 8;
            extractCasavaArgs = true;
        }
        else if (!strncmp("-", arg, 1))
        {
            extractCasavaArgs = false;
        }
        // else we're parsing a parameter argument
        if (!extractCasavaArgs)
        {
            ret.push_back(arg);
        }
        else
        {
            casavaArgv += " '";
            casavaArgv += arg;
            casavaArgv += "'";
        }
    }

    return ret;
}

void AlignOptions::parseStatsImageFormat()
{
    if(statsImageFormatString == "gif")
    {
        statsImageFormat = reports::AlignmentReportGenerator::gif;
    }
    else if(statsImageFormatString == "none")
    {
        statsImageFormat = reports::AlignmentReportGenerator::none;
    }
    else
    {
        const format message = format("\n   *** The 'stats-image-format' value is invalid %s ***\n") % statsImageFormatString;
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
}


/**
 * \brief remembers the original argv array and hands over to the base implementation
 */
common::Options::Action AlignOptions::parse(int argc, char *argv[])
{
    const std::vector<std::string> allOptions(argv, argv + argc);
    const std::vector<const char*> nonCasavaOptions = filterCasavaOptions(argv, argv + argc);
    this->argv.insert(this->argv.end(), nonCasavaOptions.begin(), nonCasavaOptions.end());
    common::Options::Action ret = common::Options::parse(nonCasavaOptions.size(), &nonCasavaOptions.front());
    if (RUN == ret)
    {
        ISAAC_THREAD_CERR << "argc: " << argc << " argv: " << boost::join(allOptions, " ") << std::endl;
    }
    return ret;
}

std::vector<flowcell::Layout::Format> AlignOptions::parseBaseCallsFormats()
{
    std::vector<flowcell::Layout::Format> ret;
    if (baseCallsFormatStringList.size() > baseCallsDirectoryList.size())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Too many --base-calls-format options specified. There must be at most one per --base-calls-directory. ***\n"));
    }
    baseCallsFormatStringList.insert(baseCallsFormatStringList.end(), baseCallsDirectoryList.size() - baseCallsFormatStringList.size(), baseCallsFormatStringList.back());
    ret.reserve(baseCallsFormatStringList.size());
    BOOST_FOREACH(const std::string &format, baseCallsFormatStringList)
    {
        if ("bcl" == format)
        {
            ret.push_back(flowcell::Layout::Bcl);
        }
        else if ("bcl-gz" == format)
        {
            ret.push_back(flowcell::Layout::BclGz);
        }
        else if ("fastq" == format)
        {
            ret.push_back(flowcell::Layout::Fastq);
        }
        else if ("fastq-gz" == format)
        {
            ret.push_back(flowcell::Layout::FastqGz);
        }
        else
        {
            BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** --base-calls-format '" + format + "' unrecognized. ***\n"));
        }
    }

    return ret;
}

std::vector<boost::filesystem::path> AlignOptions::parseSampleSheetPaths() const
{
    std::vector<boost::filesystem::path> sampleSheetPathList;
    for (size_t i = 0; baseCallsDirectoryList.size() > i; ++i)
    {
        const bfs::path &baseCallsDirectory = baseCallsDirectoryList.at(i);
        boost::filesystem::path sampleSheetPath =
            "none" == sampleSheetStringList.at(i) ? boost::filesystem::path()
          : "default" == sampleSheetStringList.at(i)
            ? boost::filesystem::exists(baseCallsDirectory / "SampleSheet.csv")
                ? baseCallsDirectory / "SampleSheet.csv" : boost::filesystem::path() : sampleSheetStringList.at(i);

        if (!sampleSheetPath.empty() && !boost::filesystem::exists(sampleSheetPath))
        {
            const boost::format message = boost::format("\n   *** Could not find supplied sample-sheet at %s ***\n") %
                sampleSheetPath;
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }
        sampleSheetPathList.push_back(sampleSheetPath);
    }
    return sampleSheetPathList;
}

void AlignOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help") ||  vm.count("version"))
    {
        return;
    }
    parseParallelization();
    verifyMandatoryPaths(vm);

    parseReferenceGenomes();

    std::vector<flowcell::Layout::Format> baseCallsFormatList = parseBaseCallsFormats();

    if (sampleSheetStringList.size() > baseCallsDirectoryList.size())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Too many --sample-sheet options specified. There must be at most one per --base-calls-directory. ***\n"));
    }
    sampleSheetStringList.insert(sampleSheetStringList.end(),
                                 baseCallsDirectoryList.size() - sampleSheetStringList.size(),
                                 "default");

    if (barcodeMismatchesStringList.size() > baseCallsDirectoryList.size())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Too many --barcode-mismatches options specified. There must be at most one per --base-calls-directory. ***\n"));
    }
    barcodeMismatchesStringList.insert(barcodeMismatchesStringList.end(),
                                 baseCallsDirectoryList.size() - barcodeMismatchesStringList.size(),
                                 barcodeMismatchesStringList.back());

    if (useBasesMaskList.size() > baseCallsDirectoryList.size())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Too many --use-bases-mask options specified. There must be at most one per --base-calls-directory. ***\n"));
    }
    useBasesMaskList.insert(useBasesMaskList.end(), baseCallsDirectoryList.size() - useBasesMaskList.size(), "default");

    if(tilesFilterList.size() > baseCallsDirectoryList.size())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Too many --tiles options specified. There must be at most one per --base-calls-directory. ***\n"));
    }
    tilesFilterList.insert(tilesFilterList.end(), baseCallsDirectoryList.size() - tilesFilterList.size(), "");

    std::vector<boost::filesystem::path> sampleSheetPathList = parseSampleSheetPaths();
    for (size_t i = 0; baseCallsDirectoryList.size() > i; ++i)
    {
        flowcell::Layout fc =
            (flowcell::Layout::Fastq == baseCallsFormatList.at(i) ||
             flowcell::Layout::FastqGz == baseCallsFormatList.at(i)) ?
            alignOptions::FastqFlowcell::createFilteredFlowcell(
                tilesFilterList.at(i),
                baseCallsDirectoryList.at(i),
                baseCallsFormatList.at(i),
                useBasesMaskList.at(i),
                seedDescriptor, referenceMetadataList) :
            alignOptions::BclFlowcell::createFilteredFlowcell(
                tilesFilterList.at(i),
                baseCallsDirectoryList.at(i),
                baseCallsFormatList.at(i),
                useBasesMaskList.at(i),
                seedDescriptor, referenceMetadataList);
        flowcellLayoutList.push_back(fc);
        // TODO: grouper does not understand FlowcellID. It expects <machine-name>_<run-number>.
        flowcellLayoutList.back().setIndex(flowcellLayoutList.size() - 1);
        const std::string flowcellId =
            (boost::format("%s_%d") %
                (fc.getFlowcellId().empty() ? "unknown-flowcell" : fc.getFlowcellId().c_str()) %
                flowcellLayoutList.back().getIndex()).str();
        flowcellLayoutList.back().setFlowcellId(flowcellId);

        const flowcell::SequencingAdapterMetadataList flowcellDefaultAdapters = defaultAdapters.size() > i ?
            parseDefaultAdapters(defaultAdapters.at(i)) : flowcell::SequencingAdapterMetadataList();
        flowcell::BarcodeMetadataList sampleSheet =
            demultiplexing::loadSampleSheetCsv(sampleSheetPathList.at(i),
                                               // use internal, guaranteed unique id for 'assumed'
                                               flowcellId,
                                               // use the id originally read from BaseCalls metadata for 'expected'
                                               fc.getFlowcellId(),
                                               flowcellLayoutList.back().getIndex(),
                                               flowcellLayoutList.back().getLaneIds(),
                                               referenceMetadataList, flowcellDefaultAdapters);

        const std::vector<unsigned> barcodeComponentMismatches = parseBarcodeMismatches(barcodeMismatchesStringList.at(i));
        std::for_each(sampleSheet.begin(), sampleSheet.end(),
                      boost::bind(&flowcell::BarcodeMetadata::setComponentMismatches, _1, boost::ref(barcodeComponentMismatches)));


        BOOST_FOREACH(flowcell::BarcodeMetadata barcode, sampleSheet)
        {
            // ensure the barcode indexes link back to their position in the global barcode list
            barcode.setIndex(barcodeMetadataList.size());
            barcodeMetadataList.push_back(barcode);
        }
    }

    if (0 >= firstPassSeeds)
    {
        const boost::format message = 
            boost::format("\n   *** At least one seed must be used on the first pass (--first-pass-seeds is %d) ***\n") % firstPassSeeds;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    parseExecutionTargets();
    parseMemoryControl();
    parseGapScoring();
    parseDodgyAlignmentScore();
    parseTemplateLength();
    parseGapRealignment();
    parseStatsImageFormat();
    parseKeepUnaligned();
}




} //namespace option
} // namespace isaac
