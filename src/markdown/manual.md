# Contents
[TOC]

# Hardware requirements

## RAM

When running from Bcl data, for human genome analyses, it is recommended to let iSAAC use at least 40 GB of RAM on a 24-threaded 
system. See [tweaks](#tweaks) section for ways to run iSAAC on limited hardware.

## IO

As a ballpark figure, if there is Y GBs of BCL data, then iSAAC roughly does the following:

* Reads 2Y GBs of BCL files
* Reads 50 GB of sorted reference (for human)
* Writes 4Y GBs of Temporary data  
* Reads 4Y GBs of Temporary Data
* Writes Y GBs of BAM Data

As a rule of thumb, given a reasonably high end modern CPU and enough memory for a lane of BCL files plus the reference, 
then the scratch storage should be able to do over 200 MB/s to avoid IO dominating the processing time and preferable over 500 MB/s.


# Sorted Reference
iSAAC requires a pre-processed reference to do the alignment. The pre-processing extracts all possible 32-mers from the 
reference genome and stores them in the format that is easily accessible to [isaac-align](#isaac-align) 
along with the metadata. The metadata file also keeps the absolute path to the original .fa file and isaac-align uses this file. 
It is important to ensure that this file is available in its original location at the time isaac-align is being run.

As the metadata uses absolute paths to reference files, manually copying or moving the sorted refernce is not recommended. 
Instead, using the [isaac-pack-reference](#isaac-pack-reference)/[isaac-unpack-reference](#isaac-unpack-reference) tool pair is advised.

In order to prepare a reference from an .fa file, use [isaac-sort-reference](#isaac-sort-reference).

# Examples

**Analyze all data from a bcl run**

    $ isaac-align -r /path/to/sorted-reference.xml -b <Run_Folder>/Data/Intensities/BaseCalls -m 40

**Analyze subset of lanes from a bcl run**

    $ isaac-align -r /path/to/sorted-reference.xml -b <Run_Folder>/Data/Intensities/BaseCalls -m 38 --tiles s_[1234]_

**Analyze paired read data from fastq.gz files**

> **NOTE:** the fastq files need to be named or sym-linked in a special way so that iSAAC can recognize them.

    $ ls Fastq/
    lane1_read1.fastq.gz  lane1_read2.fastq.gz  lane2_read1.fastq.gz  lane2_read2.fastq.gz

    $ isaac-align -r /path/to/sorted-reference.xml -b Fastq -m 40 --base-calls-format fastq-gz

**Analyze single-ended data from fastq.gz files**

> **NOTE:** that the fastq files need to be named or sym-linked in a special way so that iSAAC can recognize them.

    $ ls Fastq/
    lane1_read1.fastq.gz  lane2_read1.fastq.gz
    $ isaac-align -r /path/to/sorted-reference.xml -b Fastq -m 40 --base-calls-format fastq-gz

**Analyze data from bam file**

    $ isaac-align -r /path/to/sorted-reference.xml -b /path/to/my.bam -m 40 --base-calls-format bam

# Output folder structure

    Aligned
    |-- Projects (output data files)
    |   |-- <project name>
    |   |   |-- <sample name>
    |   |   |   |-- Casava (subset of CASAVA variant calling results data)
    |   |   |   |-- sorted.bam (bam file for the sample. Contains data for the project/sample from all flowcells)
    |   |   |   `-- sorted.bam.bai
    |   |   |-- ...
    |   `-- ...
    |-- Reports (navigable statistics pages)
    |   |-- gif
    |   |   |-- <flowcell id>
    |   |   |   `-- all
    |   |   |       `-- all
    |   |   |           `-- all
    |   |   |               |-- <per-tile statistic plot images>
    |   |   `-- ...
    |   `-- html
    |       `-- index.html (root html for the analysis reports)
    `-- Stats
        |-- BuildStats.xml (chromosome-level duplicate and coverage statistics)
        |-- DemultiplexingStats.xml (information about the barcode hits)
        `-- MatchSelectorStats.xml (tile-level yield, pair and alignment quality statistics)

# Tweaks

## Turning off data clipping

iSAAC has the following data clipping mechanisms enabled by default: [--base-quality-cutoff, --clip-semialigned and 
--clip-overlapping](#isaac-align). All have shown to improve the consistency of variant calling, however they are not suitable for 
scenarios targeting evaluation of the quality of sequencing, library preparation and such. Please see isaac-align 
command line reference manual for details.

## Reducing RAM requirement

iSAAC human genome alignment can be run in as little as 32 GB RAM by limiting the amount of compute threads with 
[-j](#isaac-align) parameter.

## Analyzing human genomes on low RAM System

As the human reference requires about 46 gigabytes to stay in RAM, the 48 gigabyte (or smaller) systems are not able to 
keep it in the cache. In this situation it pays to minimize the number of reference scanning passess the MatchFinder need to perform.

The number of tiles processed at a time is usually capped by the --temp-parallel-save parameter, which is set to 16 by 
default. This limit prevents the match files to be excessively fragmented on the systems where temporary data is stored 
on a local hard drive. If the temporary data is stored on a distributed network storage such as Isilon, the fragmentation 
is usually not an issue. In this case --temp-parallel-save 64 will ensure that entire lane (64 is the currently known 
maximum number of tiles per lane produced by an instrument) is loaded into memory for MatchFinder.

On the systems with more than 64 gigabytes of physical RAM, it makes sense to keep the defaults as in this case the 
reference has enough room to stay in the IO cache.

## Reducing Linux swappiness

iSAAC is designed to take the full advantage of the hardware resources available on the processing node. On systems with 
default Linux configuration, this causes the operating system to swap pages out to make room for the IO cache, when the 
iSAAC gets to a memory/IO intenstive stages. Since iSAAC operates on data volumes that normally exceed the amount of RAM 
of an average system, there is little or no benefit from IO-caching the data iSAAC reads and writes. As result, the data 
processing takes longer than it would if the system did not prioritize IO cache over the code and data that is already 
in RAM.

Reducing the value **_/proc/sys/vm/swappiness_** to 10 or lower from default 60 solves the problem.

# Alignment quality scoring

**Probability of a Correct Read**

This probability is the product of the probability for each base to be correct, as determined by the quality score of 
the base and the alignment of the base against the reference. If Q[i] is the Phred quality score of a base at position i, 
we have:

    pBaseError[i] = 10^(-Q[i]/10)
    pBaseCorrect[i] = 1 - pBaseError[i]
    pBaseMatch[i] = pBaseCorrect[i]
    pBaseMismatch[i] = pBaseError[i]/3
    pReadCorrect = product(i=0..readLength-1, pBase[i]), where pBase[i]is
        pBaseMatch[i] if the base at position i matches the reference, or is part of an indel
        pBaseMismatch[i] otherwise

Note: using pBaseMatch[i] for indels is an arbitrary choice (doing otherwise would require a model for indels)

**Alignment Quality of a Single Read**

The alignment quality depends on the intrinsic quality of the alignment (inferred from pReadCorrect above), but also on 
the specificity of the alignment (i.e. the probability that the read aligns somewhere else). This is inferred from two 
quantities:

    pNeighbourhood = sum of pReadCorrect for all alignments in the neighbourhood (all other identified alignment positions)
    rogCorrection = 2*GenomeLength/(4^ReadLength)

The rogCorrection is the "rest-of-genome" correction that gives an indication of the probability of having a random read 
aligning to the reference. This value tends to (and should) be very small and allows differentiating between the quality 
of reads with unique alignments (in which case pNeighbourhood == 0).

**iSAAC Accumulates Actual Probablilities in pNeighborhood**

Note that in CASAVA (unlike iSAAC), pNeighborhood is normalized to take into account all the estimated values of 
pReadCorrect for all the seeds that matched the reference with up-to two mismatches. These estimated values are 
pessimistic in the sense that they assume that extending the alignment does not introduce any additional mismatches. 
The assumption is that a seed with one mismatch will lead to an alignment descriptor with exactly one match on the base 
with the worst quality (using the definition of pReadCorrect given above). Similarly, a seed with two mismatches will 
lead to an alignment descriptor with exactly two mismatches (on the bases with the two worst qualities).

    pNormalized = rogCorrection + pNeighbourhood


finally, the alignment quality is:

    alignmentQuality == -10 * log10(pNormalized/(pNormalized + pReadCorrect))

**Alignment Quality of a Pair**

This is simply the sum of the alignment quality of each fragment when there is exactly one resolved fragment. Otherwise, 
the alignment score is corrected by the total alignment score of all the resolved fragments:

alignmentQuality = -10 * log(pBestTemplateCorrect / pTotalTemplateCorrect)

where:

    pTemplateCorrect = product(pReadCorrect for all reads)
    totalRogCorrection = rogCorrection for the total length of all reads
    pTotalTemplateCorrect = totalRogCorrection + sum(pTemplateCorrect for all resolved templates)
    pBestTemplateCorrect = max(pTemplateCorrect for all resolved templates)

# Bam files

iSAAC produces a separate bam file for each project/sample.

## Unaligned pairs

Pairs where both reads are unaligned are stored depending on the argument of [--keep-unaligned](#isaac-align) command line option.

--keep-unaligned|Behavior
:---------------|:-----------------------------------------------------------------------------------------------------
discard         |Ensures unaligned pairs are not present in the bam file
front           |Places unaligned pairs in the beginning of the bam file before the first aligned pair of the first chromosome. The iSAAC-generated bam index file is specially crafted to skip those. This approach makes it easier to locate the unaligned clusters compared to the standard implementations which require reading past the last aligned pair of the last chromosome in the genome to locate the first unaligned pair. The drawback is that the standard samtools index command is unable to process such bam files. Be sure to keep the bam index files produced by iSAAC.
back            |Makes unaligned pairs appear at the end of the bam file. Although this makes it somewhat difficult to extract unaligned data, this is the option that produced bam file that is compatible with samtools index command

## Singleton/Shadow pairs

Singleton/shadow pairs refer to pairs in which aligner was unable to decide on the alignment of one of the ends (shadow). 
In this case, the shadows are assigned the position of the end that does align (singleton). The shadows are stored in the 
bam file, immediately after their singleton.

## Read names

    read-name     = flowcell-id "_" flowcell-idx ":" lane-number ":" tile-number ":" cluster-id ":0"
    
    flowcell-id   = ;flowcell identifier from BaseCalls/config.xml. "unknown-flowcell" if the identifier cannot be 
                    ;determined from config.xml file
    flowcell-idx  = ;unique 0-based index of the flowcell within the analysis
    lane-number   = ;Lane number 1-8
    tile-number   = ;Unpadded tile number
    cluster-id    = ;Unpadded 0-based cluster id in the order in which the clusters appear in the bcl tile. 

## Bam flags usage
Bit  |Description                                            |iSAAC notes
:----|:------------------------------------------------------|:----------------------------------------------------------
0x001|template having multiple segments in sequencing        |
0x002|each segment properly aligned according to the aligner |Pair matches dominant template orientation. Single-ended templates don't have this flag set.
0x004|segment unmapped                                       |
0x008|next segment in the template unmapped                  |
0x010|SEQ being reverse complemented                         |
0x020|SEQ of the next segment in the template being reversed |
0x040|the first segment in the template                      |Read 1 (not set for single-ended data)
0x080|the last segment in the template                       |Read 2
0x100|secondary alignment                                    |iSAAC does not produce secondary alignments
0x200|not passing quality controls                           |PF flag from RTA
0x400|PCR or optical duplicate                               |If [--keep-duplicates](#isaac-align) is turned off, duplicates are excluded from the bam file. If [--mark-duplicates](#isaac-align) is turned off, duplicates are not marked in the bam file.

## Extended tags

iSAAC generates following tags in the output bam files. The list of tags stored can be controlled by --bam-exclude-tags command-line argument.

Tag|iSAAC meaning
---|:--------------------------------------------------------------------------------
AS |Pair alignment score
BC |Barcode string.
NM |Edit distance (mismatches and gaps) including the soft-clipped parts of the read
OC |Original CIGAR for the realigned reads. See --realign-gaps.
RG |iSAAC read groups correspond to flowcell/lane/barcode
SM |Single read alignment score
ZX |Cluster X pixel coordinate on the tile times 100 (disabled by default)
ZY |Cluster Y pixel coordinate on the tile times 100 (disabled by default) 

## MAPQ

Bam MAPQ for pairs that match dominant template orientation is min(max(SM, AS), 60). For reads that are not members of a 
pair machting the dominant template orientation, the MAPQ is min(SM, 60). The MAPQ might be downgraded to 0 or set to be 
unknown (255) for alignments that don't have enough evidence to be correctly scored. This behavior depends on the 
[--dodgy-alignment-score](#isaac-align) argument.

# Toolkit Reference

## isaac-align

**Usage**

    isaac-align -r <reference> -b <base calls> -m <memory limit> [optional arguments]

**Options**

    --allow-empty-flowcells arg (=0)             Avoid failure when some of the --base-calls contain no data
    --avoid-smith-waterman arg (=0)              When set, heuristics applied to avoid executing costly smith-waterman 
                                                 on sequences that are unlikely to produce gaps
    --bam-exclude-tags arg (=ZX,ZY)              Comma-separated list of regular tags to exclude from the output BAM 
                                                 files. Allowed values are: all,none,AS,BC,NM,OC,RG,SM,ZX,ZY
    --bam-gzip-level arg (=1)                    Gzip level to use for BAM
    --bam-header-tag arg                         Additional bam entries that are copied into the header of each 
                                                 produced bam file. Use '' to represent tab separators.
    --bam-pessimistic-mapq arg (=0)              When set, the MAPQ is computed as MAPQ:=min(60, min(SM, AS)), 
                                                 otherwise MAPQ:=min(60, max(SM, AS))
    --barcode-mismatches arg (=1)                Multiple entries allowed. Each entry is applied to the corresponding 
                                                 base-calls. Last entry applies to all the bases-calls-directory that 
                                                 do not have barcode-mismatches specified. Last component mismatch 
                                                 value applies to all subsequent barcode components should there be 
                                                 more than one. Examples: 
                                                   - 1:0             : allow one mismatch for the first barcode 
                                                 component and no mismatches for the subsequent components.
                                                   - 1               : allow one mismatch for every barcode component.
                                                   - 0               : no mismatches allowed in any barcode component. 
                                                 This is the default.
    -b [ --base-calls ] arg                      full path to the base calls. Multiple entries allowed. Path should 
                                                 point either to a directory or a file depending on --base-calls-format
    --base-calls-format arg                      Multiple entries allowed. Each entry is applied to the corresponding 
                                                 base-calls. Last entry is applied to all --base-calls that don't have 
                                                 --base-calls-format specified.
                                                   - bam             : --base-calls points to a Bam file. All data 
                                                 found in bam file is assumed to come from lane 1 of a single flowcell.
                                                   - bcl             : --base-calls points to RunInfo.xml file. Data is
                                                 made of uncompressed bcl files.
                                                   - bcl-gz          : --base-calls points to RunInfo.xml file. Bcl 
                                                 cycle tile files are individually compressed and named s_X_YYYY.bcl.gz
                                                   - bcl-bgzf        : --base-calls points to RunInfo.xml file. Bcl 
                                                 data is stored in cycle files that are named CCCC.bcl.bgzf
                                                   - fastq           : --base-calls points to a directory containing 
                                                 one fastq per lane/read named lane<X>_read<Y>.fastq. Use 
                                                 lane<X>_read1.fastq for single-ended data.
                                                   - fastq-gz        : --base-calls points to a directory containing 
                                                 one compressed fastq per lane/read named lane<X>_read<Y>.fastq.gz. Use
                                                 lane<X>_read1.fastq.gz for single-ended data.
    --base-quality-cutoff arg (=25)              3' end quality trimming cutoff. Value above 0 causes low quality bases
                                                 to be soft-clipped. 0 turns the trimming off.
    --bin-regex arg (=all)                       Define which bins appear in the output bam files
                                                 all                   : Include all bins in the bam and all contig 
                                                 entries in the bam header.
                                                 skip-empty             : Include only the contigs that have aligned 
                                                 data.
                                                 REGEX                 : Is treated as comma-separated list of regular 
                                                 expressions. Bam files will be filtered to contain only the bins that 
                                                 match by the name.
    --buffer-bins arg (=1)                       If set, MatchSelector will buffer bin data before writing it out. If 
                                                 not set, MatchSelector will keep an open file handle per bin and write
                                                 data into corresponding bins as it appears. This option requires extra
                                                 RAM, but improves performance on some file systems.
    --cleanup-intermediary arg (=0)              When set, iSAAC will erase intermediate input files for the stages 
                                                 that have been completed. Notice that this will prevent resumption 
                                                 from the stages that have their input files removed. --start-from Last
                                                 will still work.
    --clip-overlapping arg (=1)                  When set, the pairs that have read ends overlapping each other will 
                                                 have the lower-quality end soft-clipped.
    --clip-semialigned arg (=1)                  When set, reads have their bases soft-clipped on either sides until a 
                                                 stretch of 5 matches is found
    -c [ --cluster ] arg                         Restrict the alignment to the specified cluster Id (multiple entries 
                                                 allowed)
    --default-adapters arg                       Multiple entries allowed. Each entry is associated with the 
                                                 corresponding base-calls. Flowcells that don't have default-adapters 
                                                 provided, don't get adapters clipped in the data. 
                                                 Each entry is a comma-separated list of adapter sequences written in 
                                                 the direction of the reference. Wildcard (* character) is allowed only
                                                 on one side of the sequence. Entries with * apply only to the 
                                                 alignments on the matching strand. Entries without * apply to all 
                                                 strand alignments and are matched in the order of appearance in the 
                                                 list.
                                                 Examples:
                                                   ACGT*,*TGCA       : Will clip ACGT and all subsequent bases in the 
                                                 forward-strand alignments and mirror the behavior for the 
                                                 reverse-strand alignments.
                                                   ACGT,TGCA         : Will find the following sequences in the reads: 
                                                 ACGT, TGCA, ACGTTGCA  (but not TGCAACGT!) regardless of the alignment 
                                                 strand. Then will attempt to clip off the side of the read that is 
                                                 shorter. If both sides are roughly equal length, will clip off the 
                                                 side that has less matches.
                                                   Standard          : Standard protocol adapters. Same as 
                                                 AGATCGGAAGAGC*,*GCTCTTCCGATCT
                                                   Nextera           : Nextera standard. Same as 
                                                 CTGTCTCTTATACACATCT*,*AGATGTGTATAAGAGACAG
                                                   NexteraMp         : Nextera mate-pair. Same as 
                                                 CTGTCTCTTATACACATCT,AGATGTGTATAAGAGACAG
    --dodgy-alignment-score arg (=0)             Controls the behavior for templates where alignment score is 
                                                 impossible to assign:
                                                  - Unaligned        : marks template fragments as unaligned
                                                  - 0-254            : exact MAPQ value to be set in bam
                                                  - Unknown          : assigns value 255 for bam MAPQ. Ensures SM and 
                                                 AS are not specified in the bam
    --expected-bgzf-ratio arg (=1)               compressed = ratio * uncompressed. To avoid memory overallocation 
                                                 during the bam generation, iSAAC has to assume certain compression 
                                                 ratio. If iSAAC estimates less memory than is actually required, it 
                                                 will fail at runtime. You can check how far you are from the dangerous
                                                 zone by looking at the resident/swap memory numbers for your process 
                                                 during the bam generation. If you see too much showing as 'swap', it 
                                                 is safe to reduce the --expected-bgzf-ratio.
    --first-pass-seeds arg (=1)                  the number of seeds to use in the first pass of the match finder. Note
                                                 that this option is ignored when the --seeds=auto
    --gap-scoring arg (=bwa)                     Gapped alignment algorithm parameters:
                                                  - eland            : equivalent of 2:-1:-15:-3:-25
                                                  - bwa              : equivalent of 0:-3:-11:-4:-20
                                                  - m:mm:go:ge:me:gl : colon-delimited string of values where:
                                                      m              : match score
                                                      mm             : mismatch score
                                                      go             : gap open score
                                                      ge             : gap extend score
                                                      me             : min extend score (all gaps reaching this score 
                                                 will be treated as equal)
    --gapped-mismatches arg (=5)                 Maximum number of mismatches allowed to accept a gapped alignment.
    -h [ --help ]                                produce help message and exit
    --help-md                                    produce help message pre-formatted as a markdown file section and exit
    --ignore-missing-bcls arg (=0)               When set, missing bcl files are treated as all clusters having N bases
                                                 for the corresponding tile cycle. Otherwise, encountering a missing 
                                                 bcl file causes the analysis to fail.
    --ignore-missing-filters arg (=0)            When set, missing filter files are treated as if all clusters pass 
                                                 filter for the corresponding tile. Otherwise, encountering a missing 
                                                 filter file causes the analysis to fail.
    --ignore-neighbors arg (=0)                  When not set, MatchFinder will ignore perfect seed matches during 
                                                 single-seed pass, if the reference k-mer is known to have neighbors.
    --ignore-repeats arg (=0)                    Normally exact repeat matches prevent inexact seed matching. If this 
                                                 flag is set, inexact matches will be considered even for the seeds 
                                                 that match to repeats.
    --input-parallel-load arg (=64)              Maximum number of parallel file read operations for --base-calls
    -j [ --jobs ] arg (=32)                      Maximum number of compute threads to run in parallel
    --keep-duplicates arg (=1)                   Keep duplicate pairs in the bam file (with 0x400 flag set in all but 
                                                 the best one)
    --keep-unaligned arg (=back)                 Available options:
                                                  - discard          : discard clusters where both reads are not 
                                                 aligned
                                                  - front            : keep unaligned clusters in the front of the BAM 
                                                 file
                                                  - back             : keep unaligned clusters in the back of the BAM 
                                                 file
    --lane-number-max arg (=8)                   Maximum lane number to look for in --base-calls-directory (fastq 
                                                 only).
    --mapq-threshold arg (=0)                    Threshold used to filter the templates based on their mapping quality:
                                                 the BAM file will only contain the templates with a mapping quality 
                                                 greater than or equal to the threshold. Templates (or fragments) with 
                                                 a mapping quality of 4 or more are guaranteed to be uniquely aligned. 
                                                 Those with a mapping quality of 3 or less are either mapping to repeat
                                                 regions or have a large number of errors.
    --mark-duplicates arg (=1)                   If not set and --keep-duplicates is set, the duplicates are not 
                                                 discarded and not flagged.
    --memory-control arg (=off)                  Define the behavior in case unexpected memory allocations are 
                                                 detected: 
                                                   - warning         : Log WARNING about the allocation.
                                                   - off             : Don't monitor dynamic memory usage.
                                                   - strict          : Fail memory allocation. Intended for development
                                                 use.
    -m [ --memory-limit ] arg (=0)               Limits major memory consumption operations to a set number of 
                                                 gigabytes. 0 means no limit, however 0 is not allowed as in such case 
                                                 iSAAC will most likely consume all the memory on the system and cause 
                                                 it to crash. Default value is taken from ulimit -v.
    --neighborhood-size-threshold arg (=0)       Threshold used to decide if the number of reference 32-mers sharing 
                                                 the same prefix (16 bases) is small enough to justify the neighborhood
                                                 search. Use large enough value e.g. 10000 to enable alignment to 
                                                 positions where seeds don't match exactly.
    -o [ --output-directory ] arg (="./Aligned") Directory where the final alignment data be stored
    --output-parallel-save arg (=8)              Maximum number of parallel file write operations for 
                                                 --output-directory
    --per-tile-tls arg (=0)                      Forces template length statistics(TLS) to be recomputed for each tile.
                                                 When not set, the first tile that produces stable TLS will determine 
                                                 TLS for the rest of the tiles of the lane. Notice that as the tiles 
                                                 are not guaranteed to be processed in the same order between different
                                                 runs, some pair alignments might vary between two runs on the same 
                                                 data unless --per-tile-tls is set. It is not recommended to set 
                                                 --per-tile-tls when input data is not randomly distributed (such as 
                                                 bam) as in such cases, the shadow rescue range will be biased by the 
                                                 input data ordering.
    --pf-only arg (=1)                           When set, only the fragments passing filter (PF) are generated in the 
                                                 BAM file
    --pre-sort-bins arg (=1)                     Unset this value if you are working with references that have many 
                                                 contigs (1000+)
    --qscore-bin arg (=0)                        Toggle QScore binning, this will be applied to the data after it is 
                                                 loaded and before processing
    --qscore-bin-values arg                      Overwrite the default QScore binning values.  Default bins are 
                                                 0:0,1:1,2-9:6,10-19:15,20-24:22,25-29:27,30-34:33,35-39:37,40-63:40.  
                                                 Identity bins 1:1,2:2,3:3,4:4,5:5,6:6,7:7,8:8,9:9,10:10,11:11,12:12,13
                                                 :13,14:14,15:15,16:16,17:17,18:18,19:19,20:20,21:21,22:22,23:23,24:24,
                                                 25:25,26:26,27:27,28:28,29:29,30:30,31:31,32:32,33:33,34:34,35:35,36:3
                                                 6,37:37,38:38,39:39,40:40,41:41
    --realign-dodgy arg (=0)                     If not set, the reads without alignment score are not realigned 
                                                 against gaps found in other reads.
    --realign-gaps arg (=sample)                 For reads overlapping the gaps occurring on other reads, check if 
                                                 applying those gaps reduces mismatch count. Significantly reduces 
                                                 number of false SNPs reported around short indels.
                                                   - no              : no gap realignment
                                                   - sample          : realign against gaps found in the same sample
                                                   - project         : realign against gaps found in all samples of the
                                                 same project
                                                   - all             : realign against gaps found in all samples
    --realign-vigorously arg (=0)                If set, the realignment result will be used to search for more gaps 
                                                 and attempt another realignment, effectively extending the realignment
                                                 over multiple deletions not covered by the original alignment.
    --realigned-gaps-per-fragment arg (=1)       An estimate of how many gaps the realignment will introduce into each 
                                                 fragment.
    -r [ --reference-genome ] arg                Full path to the reference genome XML descriptor. Multiple entries 
                                                 allowed.Each entry applies to the corresponding --reference-name. The 
                                                 last --reference-genome entry may not have a corresponding 
                                                 --reference-name. In this case the default name 'default' is assumed.
    -n [ --reference-name ] arg                  Unique symbolic name of the reference. Multiple entries allowed. Each 
                                                 entry is associated with the corresponding --reference-genome and will
                                                 be matched against the 'reference' column in the sample sheet. 
                                                 Special names:
                                                   - unknown         : default reference to use with data that did not 
                                                 match any barcode.
                                                   - default         : reference to use for the data with no matching 
                                                 value in sample sheet 'reference' column.
    --repeat-threshold arg (=10)                 Threshold used to decide if matches must be discarded as too abundant 
                                                 (when the number of repeats is greater or equal to the threshold)
    -s [ --sample-sheet ] arg                    Multiple entries allowed. Each entry is applied to the corresponding 
                                                 base-calls.
                                                   - none            : process flowcell as if there is no sample sheet
                                                   - default         : use <base-calls>/SampleSheet.csv if it exists. 
                                                 This is the default behavior.
                                                   - <file path>     : use <file path> as sample sheet for the 
                                                 flowcell.
    --scatter-repeats arg (=0)                   When set, extra care will be taken to scatter pairs aligning to 
                                                 repeats across the repeat locations 
    --seed-length arg (=32)                      Length of the seed in bases. 16, 32 or 64 are allowed. Longer seeds 
                                                 reduce sensitivity on noisy data but improve repeat resolution.
    --seeds arg (=auto)                          Seed descriptors for each read, given as a comma-separated 
                                                 list-of-seeds for each read. A list-of-seeds is a colon-separated list
                                                 of offsets from the beginning of the read. 
                                                 Examples:
                                                   - auto            : automatic choice of seeds based on 
                                                 --semialigned-strategy parameter
                                                   - 0:32,0:32:64    : two seeds on the first read (at offsets 0 and 
                                                 32) and three seeds on the second read (at offsets 0, 32, and 64) and 
                                                 on subsequent reads.
                                                   - 0:32:64         : three seeds on all the reads (at offsets 0, 32 
                                                 and 64)
                                                 Note that the last list-of-seeds is repeated to all subsequent reads 
                                                 if there are more reads than there are colon-separated lists-of-seeds.
    --semialigned-gap-limit arg (=100)           The maximum length of the gap that can be introduced to minimize 
                                                 mismatches in a semialigned read. This is a separate algorithm from 
                                                 Smith-Waterman gapped alignment. use --semialigned-gap-limit 0 to 
                                                 disable this functionality.
    --shadow-scan-range arg (=-1)                -1     - scan for possible mate alignments between template min and 
                                                 max
                                                 >=0    - scan for possible mate alignments in range of template median
                                                 += shadow-scan-range
    --single-library-samples arg (=1)            If set, the duplicate detection will occur across all read pairs in 
                                                 the sample. If not set, different lanes are assumed to originate from 
                                                 different libraries and duplicate detection is not performed across 
                                                 lanes.
    --start-from arg (=Start)                    Start processing at the specified stage:
                                                   - Start            : don't resume, start from beginning
                                                   - MatchFinder      : same as Start
                                                   - MatchSelector    : skip match identification, continue with 
                                                 template selection
                                                   - AlignmentReports : regenerate alignment reports and bam
                                                   - Bam              : resume at bam generation
                                                   - Finish           : Same as Bam.
                                                   - Last             : resume from the last successful step
                                                 Note that although iSAAC attempts to perform some basic validation, 
                                                 the only safe option is 'Start' The primary purpose of the feature is 
                                                 to reduce the time required to diagnose the issues rather than be used
                                                 on a regular basis.
    --stats-image-format arg (=gif)              Format to use for images during stats generation
                                                  - gif        : produce .gif type plots
                                                  - none       : no stat generation
    --stop-at arg (=Finish)                      Stop processing after the specified stage is complete:
                                                   - Start            : perform the first stage only
                                                   - MatchFinder      : same as Start
                                                   - MatchSelector    : don't generate alignment reports and bam
                                                   - AlignmentReports : don't perform bam generation
                                                   - Bam              : finish when bam is done
                                                   - Finish           : stop at the end.
                                                   - Last             : perform up to the last successful step only
                                                 Note that although iSAAC attempts to perform some basic validation, 
                                                 the only safe option is 'Finish' The primary purpose of the feature is
                                                 to reduce the time required to diagnose the issues rather than be used
                                                 on a regular basis.
    -t [ --temp-directory ] arg (="./Temp")      Directory where the temporary files will be stored (matches, unsorted 
                                                 alignments, etc.)
    --temp-parallel-load arg (=8)                Maximum number of parallel file read operations for --temp-directory
    --temp-parallel-save arg (=64)               Maximum number of parallel file write operations for --temp-directory
    --tiles arg                                  Comma-separated list of regular expressions to select only a subset of
                                                 the tiles available in the flow-cell.
                                                 - to select all the tiles ending with '5' in all lanes: --tiles 
                                                 [0-9][0-9][0-9]5
                                                 - to select tile 2 in lane 1 and all the tiles in the other lanes: 
                                                 --tiles s_1_0002,s_[2-8]
                                                 Multiple entries allowed, each applies to the corresponding 
                                                 base-calls.
    --tls arg                                    Template-length statistics in the format 
                                                 'min:median:max:lowStdDev:highStdDev:M0:M1', where M0 and M1 are the 
                                                 numeric value of the models (0=FFp, 1=FRp, 2=RFp, 3=RRp, 4=FFm, 5=FRm,
                                                 6=RFm, 7=RRm)
    --use-bases-mask arg                         Conversion mask characters:
                                                   - Y or y          : use
                                                   - N or n          : discard
                                                   - I or i          : use for indexing
                                                 
                                                 If not given, the mask will be guessed from the B<config.xml> file in 
                                                 the base-calls directory.
                                                 
                                                 For instance, in a 2x76 indexed paired end run, the mask 
                                                 I<Y76,I6n,y75n> means:
                                                   use all 76 bases from the first end, discard the last base of the 
                                                 indexing read, and use only the first 75 bases of the second end.
    --variable-read-length arg                   Unless set, iSAAC will fail if the length of the sequence changes 
                                                 between the records of a fastq or a bam file.
    --verbosity arg (=2)                         Verbosity: FATAL(0), ERRORS(1), WARNINGS(2), INFO(3), DEBUG(4) (not 
                                                 supported yet)
    -v [ --version ]                             print program version information

## isaac-pack-reference

**Usage**

    isaac-pack-reference [options]

**Options**

    -h [ --help ]                                         Print this message
    -n [ --dry-run ]                                      Don't actually run any commands; just print them
    -v [ --version ]                                      Only print version information
    -j [ --jobs ] arg (=1)                                Maximum number of parallel operations

    -r [ --reference-genome ] arg                         Path to sorted-reference.xml 
    -o [ --output-file ] arg (./packed-reference.tar.gz)  Archive path

**Example**

    /illumina/development/iSAAC/testing/bin/isaac-pack-reference -r HumanUCSC.hg19.complete/sorted-reference.xml -j 24

## isaac-reorder-reference

**Usage**

    isaac-reorder-reference [options]

**Options**

    -b [ --bases-per-line ] arg (=70) Number of bases per line to print into .fa file.
    -h [ --help ]                     produce help message and exit
    --help-md                         produce help message pre-formatted as a markdown file section and exit
    --order arg                       Coma-separated list of contig names in the order in which they will appear in the
                                      new .fa file.
    -f [ --output-fasta ] arg         Path for the reordered fasta file.
    -x [ --output-xml ] arg           Path for the new xml file.
    -r [ --reference-genome ] arg     Full path to the reference genome XML descriptor.
    -v [ --version ]                  print program version information

## isaac-sort-reference

> **NOTE:** Available RAM could be a concern when sorting big genomes. Human genome reference sorting will require ~150 gigabytes of RAM.

**Usage**

```
isaac-sort-reference [options]
```

Time highly depends on the size of the genome. For E. coli it takes about 1 minute. For human genomes it takes about 11 
hours on a 24-threaded 2.6GHz system.

**Options**

    -g [ --genome-file ] arg                              Path to fasta file containing the reference contigs 
    -h [ --help ]                                         Print this message
    -j [ --jobs ] arg (=1)                                Maximum number of parallel operations
    -n [ --dry-run ]                                      Don't actually run any commands; just print them
    -o [ --output-directory ] arg (./iSAACIndex.<date>)   Location where the results are stored
    -q [ --quiet ]                                        Avoid excessive logging
    -p [ --no-parallel-sort ]                             Disable parallel sort when finding neighbors. Reduces RAM 
                                                          requirement by the factor of two 
    -s [ --seed-length ] arg (=32)                        Length of the k-mer. Currently 16-mer, 32-mer and 64-mer sorted 
                                                          references are supported 
    -t [ --repeat-threshold ] arg (=1000)                 Repeat cutoff after which individual kmer positions are not 
                                                          stored
    -v [ --version ]                                      Only print version information
    -w [ --mask-width ] arg (=6)                          Number of high order bits to use for splitting the 
                                                          data for parallelization 
    --dont-annotate                                       Don't search for neighbors (on by default when 
                                                          --seed-length is 16)
    --annotate                                            Force neighbor search (off by default when --seed-length 16)

**Example**

    nohup /illumina/development/iSAAC/latest/bin/isaac-sort-reference -g $(pwd)/Homo_sapiens_assembly19.fasta -j 24&

## isaac-unpack-reference

**Usage**

    isaac-unpack-reference [options]

**Options**

    -h [ --help ]                                         Print this message
    -n [ --dry-run ]                                      Don't actually run any commands; just print them
    -v [ --version ]                                      Only print version information
    -j [ --jobs ] arg (=1)                                Maximum number of parallel operations

    -w [ --mask-width ] arg (=6)                          Number of high order bits to use for splitting the 

    -i [ --input-file ] arg                               Archive path
    -t [ --repeat-threshold ] arg (=1000)                 Repeat cutoff after which individual kmer positions are not stored

