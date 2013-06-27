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
#include <boost/algorithm/string.hpp>

#include "basecalls/ConfigXml.hh"
#include "common/Exceptions.hh"
#include "demultiplexing/SampleSheetCsv.hh"
#include "oligo/Mask.hh"
#include "options/AlignOptions.hh"

#include "alignOptions/BamFlowcell.hh"
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

#define BWA_GAP_SCORING_STRING "0:-3:-11:-4:-20"
#define ELAND_GAP_SCORING_STRING "2:-1:-15:-3:-25"
static const std::vector<std::string> SUPPORTED_BAM_EXCLUDE_TAGS =
    boost::assign::list_of("AS")("BC")("NM")("OC")("RG")("SM")("ZX")("ZY");


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
    , barcodeMismatchesStringList(1, "1")
    , referenceNameList(1, "default")
    , tempDirectory("./Temp")
    , outputDirectory("./Aligned")
    , seedDescriptor("auto")//("16:0:32:64")
    , seedLength(32)
    , allowVariableReadLength(false)
    , cleanupIntermediary(false)
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
    , allowEmptyFlowcells_(false)
    , baseQualityCutoff(25) //Illumina data is normally good enough for high cutoff.
    , keepUnalignedString("back")
    , keepUnaligned(false)
    , preSortBins(true) // on genomes with decent number of contigs this speeds up bam generation by increasing memory locality
                        // of the loaded fragments. However, on metagenmoics references this causes enormous amount of entries
                        // in bin metadata data distribution
    , putUnalignedInTheBack(false)
    , realignGapsVigorously(false)
    , realignDodgyFragments(false) // true slows down pileups on DNA but seems to clear up picture significantly in RNA
    , realignedGapsPerFragment(1)
    , clipSemialigned(true) // this has to be on for some time as GATK jumps to 9000 conflict from 5000 if it is not
    , clipOverlapping(true)
    , scatterRepeats(false)
    , gappedMismatchesMax(5)
    , avoidSmithWaterman(false)
    , gapScoringString("bwa")
    , gapMatchScore(0)
    , gapMismatchScore(0)
    , gapOpenScore(0)
    , gapExtendScore(0)
    , minGapExtendScore(0)
    , semialignedGapLimit(100) //It was supposed to be 20000 or so, but GATK 1.6 and above fails with
                    //##### ERROR MESSAGE: Somehow the requested coordinate is not covered by the read. Too many deletions?
    , dodgyAlignmentScoreString("0")
    , dodgyAlignmentScore(0)
#ifdef ISAAC_THREAD_CERR_DEV_TRACE_ENABLED
    , memoryControlString("off")
#else //ISAAC_THREAD_CERR_DEV_TRACE_ENABLED
    , memoryControlString("off") //Set off by default as otherwise users get piles of WARNINGs whenever an exception is thrown.
#endif //ISAAC_THREAD_CERR_DEV_TRACE_ENABLED
    , memoryControl(common::ScoopedMallocBlock::Invalid)
    , memoryLimit(getUlimitV() / 1024 / 1024 / 1024)
    , inputLoadersMax(64) // bcl files are small, there are lots of them and at the moment they are expected to sit on a highly-parallelizable high-latency network storage
    , tempSaversMax(64)   // currently most runs using isilon as Temp storage. In this case fragmentation is not an issue
    , tempLoadersMax(8)
    , outputSaversMax(8)
    , realignGapsString("sample")
    , realignGaps(build::REALIGN_SAMPLE)
    , bamGzipLevel(boost::iostreams::gzip::best_speed)
    , expectedBgzfCompressionRatio(1)
    , singleLibrarySamples(true)
    , keepDuplicates(true)
    , markDuplicates(true)
    , binRegexString("all")
    , userTemplateLengthStatistics(-1)
    , statsImageFormatString("gif")
    , statsImageFormat(reports::AlignmentReportGenerator::gif)
    , bufferBins(true)
	, qScoreBin(false)
    , bamExcludeTags("ZX,ZY")
    , optionalFeatures(parseBamExcludeTags(bamExcludeTags))
{
    unnamedOptions_.add_options()
        ("base-calls-directory"   , bpo::value<std::vector<bfs::path> >(&baseCallsDirectoryList)->multitoken(),
            "Deprecated. Same as --base-calls.")
        ("variable-fastq-read-length"  , bpo::value<bool>(&allowVariableReadLength)->default_value(allowVariableReadLength),
            "Unless set, iSAAC will fail if the length of the sequence changes between the records of the fastq file.")
            ;

    namedOptions_.add_options()
        ("base-calls,b"   , bpo::value<std::vector<bfs::path> >(&baseCallsDirectoryList)->multitoken(),
                "full path to the base calls directory. Multiple entries allowed.")
        ("base-calls-format"        , bpo::value<std::vector<std::string> >(&baseCallsFormatStringList)->multitoken(),
                "Multiple entries allowed. Each entry is applied to the corresponding base-calls. "
                "Last entry is applied to all --base-calls that don't have --base-calls-format specified."
                "\n  - bam             : Bam file. All data found in bam file is assumed to come from lane 1 of a single "
                                        "flowcell."
                "\n  - bcl             : common bcl files, no compression."
                "\n  - bcl-gz          : bcl files are individually compressed and named s_X_YYYY.bcl.gz"
                "\n  - fastq           : One fastq per lane/read named lane<X>_read<Y>.fastq and located directly in "
                                        "the specified base-calls. Use lane<X>_read1.fastq for single-ended data."
                "\n  - fastq-gz        : One compressed fastq per lane/read named lane<X>_read<Y>.fastq.gz and "
                                        "located directly in the specified base-calls. Use lane<X>_read1.fastq.gz "
                                        "for single-ended data."
            )
        ("default-adapters"
                                    , bpo::value<std::vector<std::string> >(&defaultAdapters)->multitoken(),
                "Multiple entries allowed. Each entry is associated with the corresponding base-calls. "
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
                "\n  Standard          : Standard protocol adapters. Same as AGATCGGAAGAGC*,*GCTCTTCCGATCT"
                "\n  Nextera           : Nextera standard. Same as CTGTCTCTTATACACATCT*,*AGATGTGTATAAGAGACAG"
                "\n  NexteraMp         : Nextera mate-pair. Same as CTGTCTCTTATACACATCT,AGATGTGTATAAGAGACAG")
        ("sample-sheet,s"   , bpo::value<std::vector<std::string> >(&sampleSheetStringList)->multitoken(),
                "Multiple entries allowed. Each entry is applied to the corresponding base-calls."
                "\n  - none            : process flowcell as if there is no sample sheet"
                "\n  - default         : use <base-calls>/SampleSheet.csv if it exists. This is the default behavior."
                "\n  - <file path>     : use <file path> as sample sheet for the flowcell.")
        ("barcode-mismatches"   , bpo::value<std::vector<std::string> >(&barcodeMismatchesStringList)->multitoken()->default_value(barcodeMismatchesStringList, barcodeMismatchesStringList.at(0)),
                "Multiple entries allowed. Each entry is applied to the corresponding base-calls. "
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
                "\n  - sample          : realign against gaps found in the same sample"
                "\n  - project         : realign against gaps found in all samples of the same project"
                "\n  - all             : realign against gaps found in all samples")
        ("bam-gzip-level"           , bpo::value<int>(&bamGzipLevel)->default_value(bamGzipLevel),
                "Gzip level to use for BAM")
        ("bam-header-tag"           , bpo::value<std::vector<std::string> >(&bamHeaderTags)->multitoken(),
                "Additional bam entries that are copied into the header of each produced bam file. Use '\t' to represent tab separators.")
        ("expected-bgzf-ratio"           , bpo::value<double>(&expectedBgzfCompressionRatio)->default_value(expectedBgzfCompressionRatio),
                "compressed = ratio * uncompressed. To avoid memory overallocation during the bam generation, iSAAC has to assume certain compression ratio. "
                "If iSAAC estimates less memory than is actually required, it will fail at runtime. You can check how far "
                "you are from the dangerous zone by looking at the resident/swap memory numbers for your process "
                "during the bam generation. If you see too much showing as 'swap', it is safe to reduce the --expected-bgzf-ratio.")
        ("bam-exclude-tags"         , bpo::value<std::string>(&bamExcludeTags)->default_value(bamExcludeTags),
                ("Comma-separated list of regular tags to exclude from the output BAM files. Allowed values are: all,none," + boost::join(SUPPORTED_BAM_EXCLUDE_TAGS, ",")).c_str())
        ("tiles"                    , bpo::value<std::vector<std::string> >(&tilesFilterList),
                "Comma-separated list of regular expressions to select only a subset of the tiles available in the flow-cell."
                "\n- to select all the tiles ending with '5' in all lanes: --tiles [0-9][0-9][0-9]5"
                "\n- to select tile 2 in lane 1 and all the tiles in the other lanes: --tiles s_1_0002,s_[2-8]"
                "\nMultiple entries allowed, each applies to the corresponding base-calls.")
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
                "\n  - auto            : automatic choice of seeds based on --semialigned-strategy parameter"
                "\n  - 0:32,0:32:64    : two seeds on the first read (at offsets 0 and 32) "
                "and three seeds on the second read (at offsets 0, 32, and 64) and on subsequent reads."
                "\n  - 0:32:64         : three seeds on all the reads (at offsets 0, 32 and 64)"
                "\nNote that the last list-of-seeds is repeated to all subsequent reads "
                "if there are more reads than there are colon-separated lists-of-seeds."
            )
        ("seed-length"              , bpo::value<unsigned>(&seedLength)->default_value(seedLength),
            "Length of the seed in bases. 32 or 64 are allowed. Longer seeds reduce sensitivity on noisy data "
            "but improve repeat resolution. Ultimately 64-mer seeds are recommended for 150bp and longer data")
        ("first-pass-seeds"         , bpo::value<unsigned>(&firstPassSeeds)->default_value(firstPassSeeds),
                "the number of seeds to use in the first pass of the match finder. Note that this option "
                "is ignored when the --seeds=auto")
        ("reference-genome,r"       , bpo::value<std::vector<bfs::path> >(&sortedReferenceMetadataList),
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
                "Maximum number of parallel file read operations for --base-calls")
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
                "\n  - Finish           : Same as Bam."
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
        ("allow-empty-flowcells", bpo::value<bool>(&allowEmptyFlowcells_)->default_value(allowEmptyFlowcells_),
                "Avoid failure when some of the --base-calls contain no data")
        ("scatter-repeats"                  , bpo::value<bool>(&scatterRepeats)->default_value(scatterRepeats),
                "When set, extra care will be taken to scatter pairs aligning to repeats across the repeat locations ")
        ("base-quality-cutoff"     , bpo::value<unsigned>(&baseQualityCutoff)->default_value(baseQualityCutoff),
                "3' end quality trimming cutoff. Value above 0 causes low quality bases to be soft-clipped. 0 turns the trimming off.")
        ("variable-read-length"  , bpo::value<bool>(&allowVariableReadLength)->default_value(allowVariableReadLength),
                "Unless set, iSAAC will fail if the length of the sequence changes between the records of a fastq or a bam file.")
        ("cleanup-intermediary"  , bpo::value<bool>(&cleanupIntermediary)->default_value(cleanupIntermediary),
                "When set, iSAAC will erase intermediate input files for the stages that have been completed. Notice that "
                "this will prevent resumption from the stages that have their input files removed. --start-from Last will "
                "still work.")
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
        ("pre-sort-bins"            , bpo::value<bool>(&preSortBins)->default_value(preSortBins),
            "Unset this value if you are working with references that have many contigs (1000+)")
        ("semialigned-gap-limit"    , bpo::value<unsigned>(&semialignedGapLimit)->default_value(semialignedGapLimit),
                "The maximum length of the gap that can be introduced to minimize mismatches "
                "in a semialigned read. This is a separate algorithm from Smith-Waterman gapped alignment. "
                "use --semialigned-gap-limit 0 to disable this functionality.")
        ("clip-semialigned"         , bpo::value<bool>(&clipSemialigned)->default_value(clipSemialigned),
                "When set, reads have their bases soft-clipped on either sides until a stretch of 5 matches is found")
        ("clip-overlapping"         , bpo::value<bool>(&clipOverlapping)->default_value(clipOverlapping),
                "When set, the pairs that have read ends overlapping each other will have the lower-quality end soft-clipped.")
        ("gapped-mismatches"   , bpo::value<unsigned>(&gappedMismatchesMax)->default_value(gappedMismatchesMax),
                "Maximum number of mismatches allowed to accept a gapped alignment.")
        ("avoid-smith-waterman"     , bpo::value<bool>(&avoidSmithWaterman)->default_value(avoidSmithWaterman),
                "When set, heuristics applied to avoid executing costly smith-waterman on sequences that are unlikely to produce gaps")
        ("gap-scoring"   , bpo::value<std::string>(&gapScoringString)->default_value(gapScoringString),
                "Gapped alignment algorithm parameters:"
                "\n - eland            : equivalent of " ELAND_GAP_SCORING_STRING
                "\n - bwa              : equivalent of " BWA_GAP_SCORING_STRING
                "\n - m:mm:go:ge:me:gl : colon-delimited string of values where:"
                "\n     m              : match score"
                "\n     mm             : mismatch score"
                "\n     go             : gap open score"
                "\n     ge             : gap extend score"
                "\n     me             : min extend score (all gaps reaching this score will be treated as equal)")
        ("dodgy-alignment-score"   , bpo::value<std::string>(&dodgyAlignmentScoreString)->default_value(dodgyAlignmentScoreString),
                "Controls the behavior for templates where alignment score is impossible to assign:"
                "\n - Unaligned        : marks template fragments as unaligned"
                "\n - 0-254            : exact MAPQ value to be set in bam"
                "\n - Unknown          : assigns value 255 for bam MAPQ. Ensures SM and AS are not specified in the bam"
        )
        ("realign-vigorously"         , bpo::value<bool>(&realignGapsVigorously)->default_value(realignGapsVigorously),
                "If set, the realignment result will be used to search for more gaps and attempt another realignment, "
                "effectively extending the realignment over multiple deletions not covered by the original alignment.")
        ("realign-dodgy"         , bpo::value<bool>(&realignDodgyFragments)->default_value(realignDodgyFragments),
                "If not set, the reads without alignment score are not realigned against gaps found in other reads.")
        ("realigned-gaps-per-fragment"         , bpo::value<unsigned>(&realignedGapsPerFragment)->default_value(realignedGapsPerFragment),
                "An estimate of how many gaps the realignment will introduce into each fragment.")
        ("keep-duplicates" , bpo::value<bool>(&keepDuplicates)->default_value(keepDuplicates),
                "Keep duplicate pairs in the bam file (with 0x400 flag set in all but the best one)")
        ("mark-duplicates" , bpo::value<bool>(&markDuplicates)->default_value(markDuplicates),
                "If not set and --keep-duplicates is set, the duplicates are not discarded and not flagged.")
        ("single-library-samples" , bpo::value<bool>(&singleLibrarySamples)->default_value(singleLibrarySamples),
                "If set, the duplicate detection will occur across all read pairs in the sample. If not set, "
                "different lanes are assumed to originate from different libraries and duplicate detection is not "
                "performed across lanes.")
        ("bin-regex" , bpo::value<std::string>(&binRegexString)->default_value(binRegexString),
                "Define which bins appear in the output bam files"
                "\nall                   : Include all bins in the bam and all contig entries in the bam header."
                "\nskip-empty             : Include only the contigs that have aligned data."
                "\nREGEX                 : Is treated as comma-separated list of regular expressions. "
                "Bam files will be filtered to contain only the bins that match by the name.")
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
        ("cluster,c"                , bpo::value<std::vector<std::size_t> >(&clusterIdList)->multitoken(),
                "Restrict the alignment to the specified cluster Id (multiple entries allowed)")
        ("tls"                      , bpo::value<std::string>(&tlsString),
                "Template-length statistics in the format 'min:median:max:lowStdDev:highStdDev:M0:M1', "
                "where M0 and M1 are the numeric value of the models (0=FFp, 1=FRp, 2=RFp, 3=RRp, 4=FFm, 5=FRm, 6=RFm, 7=RRm)")
        ("stats-image-format", bpo::value<std::string>(&statsImageFormatString)->default_value(statsImageFormatString),
                "Format to use for images during stats generation"
                "\n - gif        : produce .gif type plots"
                "\n - none       : no stat generation"
        )
        ("buffer-bins"   , bpo::value<bool>(&bufferBins)->default_value(bufferBins),
                "If set, MatchSelector will buffer bin data before writing it out. If not set, MatchSelector will keep an open "
                "file handle per bin and write data into corresponding bins as it appears. This option requires extra RAM, but "
                "improves performance on some file systems.")
        ("qscore-bin"   , bpo::value<bool>(&qScoreBin)->default_value(qScoreBin),
        	    "Toggle QScore binning, this will be applied to the data after it is loaded and before processing")
        ("qscore-bin-values"   , bpo::value<std::string>(&qScoreBinValueString),
                "Overwrite the default QScore binning values.  Default bins are 0:0,1:1,2-9:6,10-19:15,20-24:22,"
        	    "25-29:27,30-34:33,35-39:37,40-63:40.  Identity bins 1:1,2:2,3:3,4:4,5:5,6:6,7:7,8:8,9:9,10:10,"
        	    "11:11,12:12,13:13,14:14,15:15,16:16,17:17,18:18,19:19,20:20,21:21,22:22,23:23,24:24,25:25,26:26,"
        		"27:27,28:28,29:29,30:30,31:31,32:32,33:33,34:34,35:35,36:36,37:37,38:38,39:39,40:40,41:41")
        ;
}

static unsigned stringToUint(const std::string &str)
{
    return boost::lexical_cast<unsigned, std::string>(str);
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
                   &stringToUint);

    return ret;
}

void AlignOptions::verifyMandatoryPaths(bpo::variables_map &vm)
{
    const std::vector<std::string> requiredOptions = boost::assign::list_of("reference-genome");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const format message = format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
    typedef std::pair<bfs::path *, std::string> PathOption;
    if (baseCallsDirectoryList.empty())
    {
        const format message = format("\n   *** At least one 'base-calls' is required ***\n");
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
    BOOST_FOREACH(bfs::path &baseCallsDirectory, baseCallsDirectoryList)
    {
        if (baseCallsDirectory.empty())
        {
            const format message = format("\n   *** The 'base-calls' can't be empty (use '.' for current directory) ***\n") % baseCallsDirectory;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        if(!exists(baseCallsDirectory))
        {
            const format message = format("\n   *** The 'base-calls' does not exist: %s ***\n") % baseCallsDirectory;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        baseCallsDirectory = boost::filesystem::absolute(baseCallsDirectory);
    }

    BOOST_FOREACH(bfs::path &sortedReferenceMetadata, sortedReferenceMetadataList)
    {
        if (sortedReferenceMetadata.empty())
        {
            const format message = format("\n   *** The 'reference-genome' can't be empty ***\n") % sortedReferenceMetadata;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        if(!exists(sortedReferenceMetadata))
        {
            const format message = format("\n   *** The 'reference-genome' does not exist: %s ***\n") % sortedReferenceMetadata;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        sortedReferenceMetadata = boost::filesystem::absolute(sortedReferenceMetadata);
        if(!exists(sortedReferenceMetadata))
        {
            const format message = format("\n   *** The 'reference-genome' does not exist: %s ***\n") % sortedReferenceMetadata;
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

build::GapRealignerMode AlignOptions::parseGapRealignment()
{
    if(realignGapsString == "yes" || realignGapsString == "sample")
    {
        return build::REALIGN_SAMPLE;
    }
    else if(realignGapsString == "project")
    {
        return build::REALIGN_PROJECT;
    }
    else if(realignGapsString == "all")
    {
        return build::REALIGN_ALL;
    }
    else if (realignGapsString != "no")
    {
        const format message = format("\n   *** The 'realign-gaps' value is invalid %s ***\n") % realignGapsString;
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
    return build::REALIGN_NONE;
}

void AlignOptions::parseExecutionTargets()
{
    const static std::vector<std::string> allowedStageStrings =
        boost::assign::list_of("Start")("MatchFinder")("MatchSelector")
            ("AlignmentReports")("Bam")
            ("Finish")("Last");

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
        5 == stopAtPos ? workflow::AlignWorkflow::Finish :
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
        gapScoringString = BWA_GAP_SCORING_STRING;
    }
    else if ("eland" == gapScoringString)
    {
        gapScoringString = ELAND_GAP_SCORING_STRING;
    }

    {
        std::vector<std::string> gapScores;
        using boost::algorithm::split;
        using boost::algorithm::is_any_of;

        boost::split(gapScores, gapScoringString,  is_any_of(":"));

        if (5 != gapScores.size())
        {
            const format message = format("\n   *** The 'gap-scoring' string must contain five components delimited by ':' ***\n");
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }

        gapMatchScore = boost::lexical_cast<int>(gapScores.at(0));
        if (0 > gapMatchScore)
        {
            const format message = format("\n   *** The 'gap-scoring' string must contain positive value or 0 for match score ***\n");
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        gapMismatchScore = boost::lexical_cast<int>(gapScores.at(1));
        if (0 < gapMismatchScore)
        {
            const format message = format("\n   *** The 'gap-scoring' string must contain negative value or 0 for mismatch score ***\n");
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        gapOpenScore = boost::lexical_cast<int>(gapScores.at(2));
        if (0 < gapOpenScore)
        {
            const format message = format("\n   *** The 'gap-scoring' string must contain negative value or 0 for gap open score ***\n");
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        gapExtendScore = boost::lexical_cast<int>(gapScores.at(3));
        if (0 < gapExtendScore)
        {
            const format message = format("\n   *** The 'gap-scoring' string must contain negative value or 0 for gap extend score ***\n");
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
        minGapExtendScore = boost::lexical_cast<int>(gapScores.at(4));
        if (0 < minGapExtendScore)
        {
            const format message = format("\n   *** The 'gap-scoring' string must contain negative value or 0 for gap extend score cap ***\n");
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
}

workflow::AlignWorkflow::OptionalFeatures AlignOptions::parseBamExcludeTags(std::string strBamExcludeTags)
{
    // note this must be one entry per workflow::AlignWorkflow::OptionalFeatures in the same order
    workflow::AlignWorkflow::OptionalFeatures ret = workflow::AlignWorkflow::Everything;
    if ("all" == strBamExcludeTags)
    {
        ret = workflow::AlignWorkflow::Nothing;
    }
    else if ("none" != strBamExcludeTags)
    {
        std::vector<std::string> bamExcludeTagsList;
        boost::algorithm::split(bamExcludeTagsList, strBamExcludeTags,  boost::algorithm::is_any_of(","));
        std::sort(bamExcludeTagsList.begin(), bamExcludeTagsList.end());

        BOOST_FOREACH(const std::string &excludeTag, bamExcludeTagsList)
        {
            const std::vector<std::string>::const_iterator it = std::find(SUPPORTED_BAM_EXCLUDE_TAGS.begin(), SUPPORTED_BAM_EXCLUDE_TAGS.end(), excludeTag);
            if (SUPPORTED_BAM_EXCLUDE_TAGS.end() == it)
            {
                const format message = format("\n   *** The '--bam-exclude-tags' contains unsupported tag: %s. Only following are allowed: %s ***\n") %
                    excludeTag % boost::join(SUPPORTED_BAM_EXCLUDE_TAGS, ",");
                BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
            }
            ret = workflow::AlignWorkflow::OptionalFeatures(ret & ~(1 << std::distance(SUPPORTED_BAM_EXCLUDE_TAGS.begin(), it)));
        }
    }
    return ret;
}

void AlignOptions::parseBamExcludeTags()
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
        dodgyAlignmentScore = alignment::TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNKNOWN;
    }
    else if ("Unaligned" == dodgyAlignmentScoreString)
    {
        dodgyAlignmentScore = alignment::TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED;
    }
    else
    {
    	try
    	{
    		dodgyAlignmentScore = boost::lexical_cast<short>(dodgyAlignmentScoreString);
    	}
    	catch (boost::bad_lexical_cast &e)
    	{
    		const format message = format("\n   *** The 'dodgy-alignment-score' option must be either Unknown, Unaligned or a number 0-255 (%s given) ***\n") % dodgyAlignmentScoreString;
    		BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    	}

    	if (dodgyAlignmentScore < 0 || dodgyAlignmentScore > 254)
    	{
    		const format message = format("\n   *** The 'dodgy-alignment-score' option must be either Unknown, Unaligned or a number 0-255 (%s given) ***\n") % dodgyAlignmentScoreString;
    		BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    	}
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
    if (referenceNameList.size() > sortedReferenceMetadataList.size())
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

    if (referenceNameList.size() < sortedReferenceMetadataList.size())
    {
        if (referenceNameList.end() == std::find(referenceNameList.begin(), referenceNameList.end(), "default"))
        {
            referenceNameList.push_back("default");
        }
    }

    if (referenceNameList.size() < sortedReferenceMetadataList.size())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Not enough --reference-name options specified. There must be one for each --reference-genome. ***\n"));
    }

    BOOST_FOREACH(const boost::filesystem::path &sortedReferenceMetadata, sortedReferenceMetadataList)
    {
        const unsigned referenceIndex = &sortedReferenceMetadata - &sortedReferenceMetadataList.front();
        referenceMetadataList.push_back(
            reference::ReferenceMetadata(referenceNameList.at(referenceIndex), sortedReferenceMetadata,
                                         referenceMetadataList.size()));
    }
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
    ISAAC_THREAD_CERR << "Version: " << iSAAC_VERSION_FULL << std::endl;
    ISAAC_THREAD_CERR << "argc: " << argc << " argv: " << boost::join(allOptions, " ") << std::endl;

    this->argv.insert(this->argv.end(), argv, argv + argc);
    common::Options::Action ret = common::Options::parse(argc, argv);
    return ret;
}

std::vector<flowcell::Layout::Format> AlignOptions::parseBaseCallsFormats()
{
    std::vector<flowcell::Layout::Format> ret;
    if (baseCallsFormatStringList.size() > baseCallsDirectoryList.size())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Too many --base-calls-format options specified. There must be at most one per --base-calls. ***\n"));
    }
    baseCallsFormatStringList.insert(baseCallsFormatStringList.end(), baseCallsDirectoryList.size() - baseCallsFormatStringList.size(), baseCallsFormatStringList.back());
    ret.reserve(baseCallsFormatStringList.size());
    BOOST_FOREACH(const std::string &format, baseCallsFormatStringList)
    {
        if ("bam" == format)
        {
            ret.push_back(flowcell::Layout::Bam);
        }
        else if ("bcl" == format)
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
    for (std::size_t i = 0; baseCallsDirectoryList.size() > i; ++i)
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

flowcell::BarcodeMetadataList::const_iterator validateBarcodeGroupReferenceMapping(
    const flowcell::BarcodeMetadataList::const_iterator groupBegin,
    const flowcell::BarcodeMetadataList::const_iterator groupEnd)
{
    unsigned reference = -1;
    // check that all barcodes are aligned against the same reference
    for(flowcell::BarcodeMetadataList::const_iterator it = groupBegin; groupEnd != it; ++it)
    {
        const flowcell::BarcodeMetadata &barcode = *it;
        if (!barcode.isUnmappedReference())
        {
            if (-1U == reference)
            {
                reference = barcode.getReferenceIndex();
            }
            else
            {
                if (reference != barcode.getReferenceIndex())
                {
                    return it;
                }
            }
        }
    }
    return groupEnd;
}

inline bool orderByBarcodeProject(
    const flowcell::BarcodeMetadata &left,
    const flowcell::BarcodeMetadata &right)
{
    return left.getProject() < right.getProject();
}

inline bool orderByBarcodeSample(
    const flowcell::BarcodeMetadata &left,
    const flowcell::BarcodeMetadata &right)
{
    return left.getSampleName() < right.getSampleName();
}

void validateSampleSheets(build::GapRealignerMode realignGaps, flowcell::BarcodeMetadataList barcodes)
{
    if (build::REALIGN_NONE == realignGaps)
    {
        // no validation is required
        return;
    }
    else if (build::REALIGN_ALL == realignGaps)
    {
        const flowcell::BarcodeMetadataList::const_iterator invalidReferenceBcIt =
            validateBarcodeGroupReferenceMapping(barcodes.begin(), barcodes.end());
        if (barcodes.end() != invalidReferenceBcIt)
        {
            const boost::format message = boost::format("\n   *** "
                "Barcode reference mapping is incompatible with the rest of the barcodes: %s. "
                "All barcodes must map to the same reference when '--realign-gaps all' is used ***\n") %
                *invalidReferenceBcIt;
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }
    }
    else if (build::REALIGN_PROJECT == realignGaps)
    {
        std::sort(barcodes.begin(), barcodes.end(), orderByBarcodeProject);
        flowcell::BarcodeMetadataList::const_iterator currentProjectBcIt = barcodes.begin();
        while(barcodes.end() != currentProjectBcIt)
        {
            const flowcell::BarcodeMetadataList::const_iterator nextProjectBcIt =
                std::upper_bound<flowcell::BarcodeMetadataList::const_iterator>(currentProjectBcIt + 1, barcodes.end(),
                                                                                *currentProjectBcIt, orderByBarcodeProject);

            const flowcell::BarcodeMetadataList::const_iterator invalidReferenceBcIt =
                validateBarcodeGroupReferenceMapping(currentProjectBcIt, nextProjectBcIt);

            if (nextProjectBcIt != invalidReferenceBcIt)
            {
                const boost::format message = boost::format("\n   *** "
                    "Barcode reference mapping is incompatible with the rest of the barcodes for project %s: %s. "
                    "All barcodes of the project must map to the same reference when '--realign-gaps project' is used ***\n") %
                    currentProjectBcIt->getProject() % *invalidReferenceBcIt;
                BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
            }
            currentProjectBcIt = nextProjectBcIt;
        }
    }
    else if (build::REALIGN_SAMPLE == realignGaps)
    {
        std::sort(barcodes.begin(), barcodes.end(), orderByBarcodeSample);
        flowcell::BarcodeMetadataList::const_iterator currentSampleBcIt = barcodes.begin();
        while(barcodes.end() != currentSampleBcIt)
        {
            const flowcell::BarcodeMetadataList::const_iterator nextSampleBcIt =
                std::upper_bound<flowcell::BarcodeMetadataList::const_iterator>(currentSampleBcIt + 1, barcodes.end(),
                                                                                *currentSampleBcIt, orderByBarcodeSample);

            const flowcell::BarcodeMetadataList::const_iterator invalidReferenceBcIt =
                validateBarcodeGroupReferenceMapping(currentSampleBcIt, nextSampleBcIt);

            if (nextSampleBcIt != invalidReferenceBcIt)
            {
                const boost::format message = boost::format("\n   *** "
                    "Barcode reference mapping is incompatible with the rest of the barcodes for sample %s: %s. "
                    "All barcodes of the project must map to the same reference when '--realign-gaps project' is used ***\n") %
                    currentSampleBcIt->getSampleName() % *invalidReferenceBcIt;
                BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
            }
            currentSampleBcIt = nextSampleBcIt;
        }
    }
}

inline void unescapeSlashT(std::string &str)
{
    for (std::size_t pos = str.find("\\t"); std::string::npos != pos; pos = str.find("\\t", pos))
    {
        str.replace(pos, 2, "\t");
    }
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
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Too many --sample-sheet options specified. There must be at most one per --base-calls. ***\n"));
    }
    sampleSheetStringList.insert(sampleSheetStringList.end(),
                                 baseCallsDirectoryList.size() - sampleSheetStringList.size(),
                                 "default");

    if (barcodeMismatchesStringList.size() > baseCallsDirectoryList.size())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Too many --barcode-mismatches options specified. There must be at most one per --base-calls. ***\n"));
    }
    barcodeMismatchesStringList.insert(barcodeMismatchesStringList.end(),
                                 baseCallsDirectoryList.size() - barcodeMismatchesStringList.size(),
                                 barcodeMismatchesStringList.back());

    if (useBasesMaskList.size() > baseCallsDirectoryList.size())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Too many --use-bases-mask options specified. There must be at most one per --base-calls. ***\n"));
    }
    useBasesMaskList.insert(useBasesMaskList.end(), baseCallsDirectoryList.size() - useBasesMaskList.size(), "default");

    if(tilesFilterList.size() > baseCallsDirectoryList.size())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Too many --tiles options specified. There must be at most one per --base-calls. ***\n"));
    }
    tilesFilterList.insert(tilesFilterList.end(), baseCallsDirectoryList.size() - tilesFilterList.size(), "");

    if (semialignedGapLimit)
    {
        if ("auto" == seedDescriptor)
        {
            firstPassSeeds = 2;
        }
    }

    if (16 != seedLength && 32 != seedLength && 64 != seedLength)
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** --seed-length other than 16, 32 or 64 is not supported. ***\n"));
    }

    std::vector<boost::filesystem::path> sampleSheetPathList = parseSampleSheetPaths();
    for (std::size_t i = 0; baseCallsDirectoryList.size() > i; ++i)
    {
        flowcell::Layout fc =
            (flowcell::Layout::Fastq == baseCallsFormatList.at(i) ||
             flowcell::Layout::FastqGz == baseCallsFormatList.at(i)) ?
            alignOptions::FastqFlowcell::createFilteredFlowcell(
                semialignedGapLimit,
                tilesFilterList.at(i),
                baseCallsDirectoryList.at(i),
                baseCallsFormatList.at(i),
                useBasesMaskList.at(i),
                allowVariableReadLength,
                seedDescriptor, seedLength, referenceMetadataList,
                firstPassSeeds) :
            flowcell::Layout::Bam == baseCallsFormatList.at(i) ?
            alignOptions::BamFlowcell::createFilteredFlowcell(
                semialignedGapLimit,
                tilesFilterList.at(i),
                baseCallsDirectoryList.at(i),
                baseCallsFormatList.at(i),
                useBasesMaskList.at(i),
                allowVariableReadLength,
                seedDescriptor, seedLength, referenceMetadataList,
                firstPassSeeds) :
            alignOptions::BclFlowcell::createFilteredFlowcell(
                semialignedGapLimit,
                tilesFilterList.at(i),
                baseCallsDirectoryList.at(i),
                baseCallsFormatList.at(i),
                useBasesMaskList.at(i),
                seedDescriptor, seedLength, referenceMetadataList,
                firstPassSeeds);
        if (!fc.getReadMetadataList().empty())
        {
            if (fc.getLaneIds().empty())
            {
                if (!allowEmptyFlowcells_)
                {
                    const boost::format message = boost::format("\n   *** Could not find any lanes matching the '%s' in: %s. Use --allow-empty-flowcell to avoid the failure. ***\n") %
                        tilesFilterList.at(i) % baseCallsDirectoryList.at(i);
                    BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
                }
            }
            else
            {
                const std::string oriFlowcellId = fc.getFlowcellId().empty() ? "unknown-flowcell" : fc.getFlowcellId().c_str();
                std::string uniqueFlowcellId = oriFlowcellId;
                unsigned conflictingIdCount = 0;
                while (flowcellLayoutList.end() !=
                    std::find_if(flowcellLayoutList.begin(), flowcellLayoutList.end(),
                                 boost::bind(&flowcell::Layout::getFlowcellId, _1) == uniqueFlowcellId))
                {
                    ++conflictingIdCount;
                    uniqueFlowcellId = oriFlowcellId + (boost::format("-%d") % conflictingIdCount).str();
                    ISAAC_THREAD_CERR << uniqueFlowcellId << std::endl;
                }
                if (oriFlowcellId != uniqueFlowcellId)
                {
                    ISAAC_THREAD_CERR << "WARNING: renamed flowcell id " << oriFlowcellId << " into " << uniqueFlowcellId << " to avoid duplication" << std::endl;
                }

                fc.setFlowcellId(uniqueFlowcellId);

                flowcellLayoutList.push_back(fc);
                // TODO: grouper does not understand FlowcellID. It expects <machine-name>_<run-number>.
                flowcellLayoutList.back().setIndex(flowcellLayoutList.size() - 1);

                const flowcell::SequencingAdapterMetadataList flowcellDefaultAdapters = defaultAdapters.size() > i ?
                    parseDefaultAdapters(defaultAdapters.at(i)) : flowcell::SequencingAdapterMetadataList();
                flowcell::BarcodeMetadataList sampleSheet =
                    demultiplexing::loadSampleSheetCsv(sampleSheetPathList.at(i),
                                                       // use internal, guaranteed unique id for 'assumed'
                                                       uniqueFlowcellId,
                                                       // use the id originally read from BaseCalls metadata for 'expected'
                                                       fc.getFlowcellId(),
                                                       flowcellLayoutList.back().getBarcodeLength(),
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
        }
        else
        {
            if (allowEmptyFlowcells_)
            {
                ISAAC_THREAD_CERR << "WARNING: Flowcell " << fc << " has no data" << std::endl;
            }
            else
            {
                const boost::format message =
                    boost::format("\n   *** %s has no data. Use --allow-empty-flowcell to avoid the failure. ***\n") % fc;
                BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
            }
        }

        if (flowcellLayoutList.empty())
        {
            BOOST_THROW_EXCEPTION(common::InvalidOptionException("No data found to process. Please check your --base-calls."));
        }
    }

    if (0 >= firstPassSeeds)
    {
        const boost::format message = 
            boost::format("\n   *** At least one seed must be used on the first pass (--first-pass-seeds is %d) ***\n") % firstPassSeeds;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    realignGaps = parseGapRealignment();
    std::for_each(bamHeaderTags.begin(), bamHeaderTags.end(), unescapeSlashT);
    validateSampleSheets(realignGaps, barcodeMetadataList);

    parseExecutionTargets();
    parseMemoryControl();
    parseGapScoring();
    parseDodgyAlignmentScore();
    parseTemplateLength();
    parseStatsImageFormat();
    optionalFeatures = parseBamExcludeTags(bamExcludeTags);
    parseQScoreBinValues();
    parseBamExcludeTags();
}


// Set the score of the boost::array
void setScore(boost::array<char, 256> &table, const unsigned int idx, const unsigned int value)
{
	for(unsigned int base = 0; base < 4; ++base)
	{
		table[(idx << 2) | base] = (value << 2) | base;
	}
}


void AlignOptions::parseQScoreBinValues()
{
	// Check to make sure there are sane settings
	if(!qScoreBin && qScoreBinValueString.size())
		BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** QScore bin set to off but bin values supplied, remove bin values ***\n"));

	// First initialize to default bins
	const int MAXSCORE = 63;
	for(int s = 0; s < MAXSCORE+1; ++s)
	{
		setScore(fullBclQScoreTable, s,
			(0 == s) ? 0 :
			(1 == s) ? 1 :
			(10 > s) ? 6 :
			(20 > s) ? 15:
			(25 > s) ? 22:
			(30 > s) ? 27:
			(35 > s) ? 33:
			(40 > s) ? 37:
			40);
	}

	// Now that full table is created just return if there is nothing to do
	if(qScoreBinValueString.empty()) return;

	// Validate user input
	const std::string validQScoreRegexStr =  "^(?:(?:(\\d+)-(\\d+)):(\\d+),|(\\d+):(\\d+),)*(?:(?:(\\d+)-(\\d+)):(\\d+)|(\\d+):(\\d+))$";
	const boost::regex validQScoreRegex(validQScoreRegexStr);
	boost::smatch what;
	if(!boost::regex_match(qScoreBinValueString, what, validQScoreRegex))
	{
		BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Error validating qscore bin ranges " + qScoreBinValueString + " ***\n"));
	}

	// Parse user input
	const std::string binStr = "(?:(\\d+)-(\\d+)|(\\d+)):(\\d+)";
	// For each match in qScoreBinValueString
	BOOST_FOREACH(const boost::match_results<std::string::const_iterator> &what,
		boost::make_iterator_range(boost::sregex_iterator(qScoreBinValueString.begin(),qScoreBinValueString.end(),boost::regex(binStr)),boost::sregex_iterator()))
	{
		if(what.size() != 5)
			BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Error with qScore regex parsing " + what.str() + ", capture size incorrect ***\n"));

		// Get ranges and bin score
		const int rangea = boost::lexical_cast<int>( (what[1].second != what[1].first) ?
			std::string(what[1].first, what[1].second) :
			std::string(what[3].first, what[3].second) );
		const int rangeb =  (what[2].second != what[2].first) ?
			boost::lexical_cast<int>(std::string(what[2].first, what[2].second)) :
			rangea;
		const int binScore = boost::lexical_cast<int>(std::string(what[4].first, what[4].second));

		// Make sure range and bin score are in acceptable range
		if(rangea > MAXSCORE || rangea < 0)
			BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Error with qScore bin start range " + boost::lexical_cast<std::string>(rangea) + ", must be 0-63 ***\n"));
		if(rangeb > MAXSCORE || rangea < 0)
			BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Error with qScore bin end range " + boost::lexical_cast<std::string>(rangeb) + ", must be 0-63 ***\n"));
		if(binScore > MAXSCORE || rangea < 0)
			BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Error with qScore bin score " + boost::lexical_cast<std::string>(binScore) + ", must be 0-63 ***\n"));

		// Set range in table
		for(int i = rangea; i <= rangeb; ++i)
		{
			setScore(fullBclQScoreTable, i, binScore);
		}
	}
}


} //namespace option
} // namespace isaac
