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
 ** \file Build.cpp
 **
 ** Reorders aligned data and stores results in bam file.
 **
 ** \author Roman Petrovski
 **/

#include "common/config.h"

#ifdef HAVE_NUMA
#include <numa.h>
#include <numaif.h>
#endif //HAVE_NUMA
 
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include "bam/Bam.hh"
#include "bam/BamIndexer.hh"
#include "bam/BgzfCompressor.hh"
#include "build/Build.hh"
#include "common/Debug.hh"
#include "common/FileSystem.hh"
#include "common/Threads.hpp"
#include "io/Fragment.hh"
#include "reference/ContigLoader.hh"

#include "BuildStatsXml.hh"
#include "SortedReferenceXmlBamHeaderAdapter.hh"

namespace isaac
{
namespace build
{

/**
 * \return Returns the total memory in bytes required to load the bin data and indexes
 */
static unsigned long getBinTotalSize(const alignment::BinMetadata & binMetadata)
{
    return
        binMetadata.getDataSize() +
        binMetadata.getFIdxElements() * sizeof(io::FStrandFragmentIndex) +
        binMetadata.getRIdxElements() * sizeof(io::RStrandOrShadowFragmentIndex) +
        binMetadata.getSeIdxElements() * sizeof(io::SeFragmentIndex);
}

unsigned long Build::estimateBinCompressedDataRequirements(
    const alignment::BinMetadata & binMetadata,
    const unsigned outputFileIndex) const
{
    // TODO: puth the real number in here.
    static const unsigned long EMPTY_BGZF_BLOCK_SIZE = 1234UL;
    if (!binMetadata.getTotalElements())
    {
        return EMPTY_BGZF_BLOCK_SIZE;
    }

    unsigned long thisOutputFileBarcodeElements = 0;
    // accumulate size required to store all barcodes that map to the same output file
    BOOST_FOREACH(const flowcell::BarcodeMetadata& barcode, barcodeMetadataList_)
    {
        const unsigned barcodeOutputFileIndex = barcodeBamMapping_.getFileIndex(barcode);
        if (outputFileIndex == barcodeOutputFileIndex)
        {
            const unsigned barcodeIndex = barcode.getIndex();
            ISAAC_ASSERT_MSG(0 != binMetadata.getTotalElements() || 0 == binMetadata.getBarcodeElements(barcodeIndex), "Can't have empty bin with non-empty bin barcode");

            thisOutputFileBarcodeElements += binMetadata.getBarcodeElements(barcodeIndex);
        }
    }

    // assume all data will take the same fraction or less than the number derived from demultiplexed fragments.
    return EMPTY_BGZF_BLOCK_SIZE +
        ((getBinTotalSize(binMetadata) * thisOutputFileBarcodeElements +
            binMetadata.getTotalElements() - 1) / binMetadata.getTotalElements()) * expectedBgzfCompressionRatio_;
}

/**
 * \brief Produces mapping so that all barcodes having the same sample name go into the same output file.
 *
 * \return Returns the pair of a vector of that maps a barcode index to a unique output file index in the
 *         second vector so that two barcodes that are supposed to go into the same file will end up
 *         having the same mapping
 */
BarcodeBamMapping mapBarcodesToFiles(
    const boost::filesystem::path &outputDirectory,
    const flowcell::BarcodeMetadataList &barcodeMetadataList)
{
    BarcodeBamMapping ret;
    unsigned outputFileIndex = 0;
    typedef std::map<boost::filesystem::path, unsigned > PathToIndex;
    std::vector<boost::filesystem::path> directories;
    PathToIndex pathToIdx;
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        directories.push_back(outputDirectory / barcode.getProject());
        directories.push_back(outputDirectory / barcode.getProject() / barcode.getSampleName());
        boost::filesystem::path bamPath(directories.back() / "sorted.bam");

        PathToIndex::iterator it = pathToIdx.find(bamPath);
        if (pathToIdx.end() == it)
        {
            it = pathToIdx.insert(std::make_pair(bamPath, outputFileIndex++)).first;
            ret.mapToNew(barcode, bamPath);
        }
        else
        {
            ret.mapToExisting(barcode, it->second);
        }
    }
    common::createDirectories(directories);

    return ret;
}

std::vector<boost::shared_ptr<boost::iostreams::filtering_ostream> > Build::createOutputFileStreams(
    const flowcell::TileMetadataList &tileMetadataList,
    const flowcell::BarcodeMetadataList &barcodeMetadataList)
{
    unsigned sinkIndexToCreate = 0;
    std::vector<boost::shared_ptr<boost::iostreams::filtering_ostream> > ret;
    ret.reserve(barcodeBamMapping_.getTotalFiles());

    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        if (sinkIndexToCreate == barcodeBamMapping_.getFileIndex(barcode))
        {
            const boost::filesystem::path &bamPath = barcodeBamMapping_.getFilePath(barcode);
            if (!barcode.isUnmappedReference())
            {
                boost::iostreams::file_sink sink(bamPath.string());
                if (!sink.is_open()) {
                    BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open output BAM file " + bamPath.string()));
                }
                ISAAC_THREAD_CERR << "Created BAM file: " << bamPath << std::endl;

                ret.push_back(boost::shared_ptr<boost::iostreams::filtering_ostream>(new boost::iostreams::filtering_ostream));
                boost::iostreams::filtering_ostream &bamStream = *ret.back();

//#define USE_OLD_INDEXER
#ifdef USE_OLD_INDEXER
                // Add BAM Indexer
                const boost::filesystem::path &baiPath = bamPath.string() + ".bai.old";
                boost::iostreams::file_sink baiSink(baiPath.string());
                if (!baiSink.is_open()) {
                    BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open output BAI file " + baiPath.string()));
                }
                ISAAC_THREAD_CERR << "Created BAI file: " << baiPath << std::endl;
                bamStream.push(isaac::bam::BamIndexer<boost::iostreams::file_sink>(baiSink));
#endif //ifdef USE_OLD_INDEXER
                bamStream.push(sink);

                {
                    boost::iostreams::filtering_ostream bgzfStream;
                    bgzfStream.push(bam::BgzfCompressor(bamGzipLevel_));
                    bgzfStream.push(bamStream);

                    bam::serializeHeader(bgzfStream,
                                         argv_,
                                         SortedReferenceXmlBamHeaderAdapter(
                                             sortedReferenceXmlList_.at(barcode.getReferenceIndex()),
                                             tileMetadataList, barcodeMetadataList,
                                             barcode.getSampleName()));
                }

                // Create BAM Indexer
                bamStream.strict_sync();
#ifdef USE_OLD_INDEXER
                // If the old indexer is activated, we can't do bamStream.tellp
                // This alternative is not very safe, but safe enough for debugging
                unsigned headerCompressedLength = boost::filesystem::file_size(bamPath);
                std::clog << "headerCompressedLength=" << headerCompressedLength << std::endl;
#else //ifdef USE_OLD_INDEXER
                unsigned headerCompressedLength = bamStream.tellp();
#endif //ifdef USE_OLD_INDEXER
                unsigned contigCount = sortedReferenceXmlList_.at(barcode.getReferenceIndex()).getContigsCount();
                bamIndexes_.push_back(new bam::BamIndex(bamPath, contigCount, headerCompressedLength));
            }
            else
            {
                // dimensions must match. Also, the way the unaligned records are implemented, requires a valid stream.
                ret.push_back(boost::shared_ptr<boost::iostreams::filtering_ostream>(new boost::iostreams::filtering_ostream));
                boost::iostreams::filtering_ostream &bamStream = *ret.back();
                bamStream.push(boost::iostreams::null_sink());
                bamIndexes_.push_back(new bam::BamIndex());
                ISAAC_THREAD_CERR << "Skipped BAM file due to unmapped barcode reference: " << bamPath << " " << barcode << std::endl;
            }
            ++sinkIndexToCreate;
        }
    }
    ISAAC_ASSERT_MSG(barcodeBamMapping_.getTotalFiles() == sinkIndexToCreate, "must create all output file sinks");

    return ret;
}

static alignment::BinMetadataList filterBins(
    const alignment::BinMetadataList& bins,
    const std::string &binRegexString)
{
    if (binRegexString.empty())
    {
        return bins;
    }
    alignment::BinMetadataList ret;
    std::string regexString(binRegexString);
    std::replace(regexString.begin(), regexString.end(), ',', '|');
    boost::regex re(regexString);
    BOOST_FOREACH(const alignment::BinMetadata &bin, bins)
    {
        if (boost::regex_search(bin.getPath().filename().string(), re))
        {
            ret.push_back(bin);
        }
    }
    if (ret.empty())
    {
        ISAAC_THREAD_CERR << "WARNING: Bam files will be empty. No bins are left after applying the following regex filter: "
            << regexString << std::endl;
    }
    return ret;
}

static void breakUpBin(
    const alignment::BinMetadata& bin,
    const unsigned long partsCount,
    alignment::BinMetadataList &ret)
{
    ISAAC_ASSERT_MSG(0 == bin.getIndex(), "At the moment only bin 0 is expected to be unaligned");
    const unsigned long newBinSize = bin.getDataSize() / partsCount;
    ISAAC_THREAD_CERR << "Breaking unaligned bin of " << bin.getDataSize()/1024/1024 << " megabytes into " <<
        partsCount << " bins of " << newBinSize /1024/1024 << " megabytes for better parallelization: " << bin << std::endl;

    for (unsigned long offset = 0; bin.getDataSize() > offset;)
    {
        alignment::BinMetadata part = bin.getChunks(offset, newBinSize);
        ISAAC_THREAD_CERR << " offset:" << offset << " " << part <<std::endl;
        offset += part.getDataSize();
        part.setIndex(ret.size());
        ret.push_back(part);
    }
}

/**
 * \brief breaks unaligned bin into about partsCount of roughly equivalent size bins
 */
static alignment::BinMetadataList breakUpUnalignedBin(
    const alignment::BinMetadataList& bins,
    const unsigned long partsCount,
    const bool keepUnaligned,
    const bool putUnalignedInTheBack)
{
    if (1 == partsCount && 1 >= bins.size())
    {
        return bins;
    }
    alignment::BinMetadataList ret;

    // put all the unaligned in the front if required
    if (keepUnaligned && !putUnalignedInTheBack)
    {
        BOOST_FOREACH(alignment::BinMetadata bin, bins)
        {
            if (bin.isUnalignedBin())
            {
                breakUpBin(bin, partsCount, ret);
            }
        }
    }

    // put all the regular ones
    BOOST_FOREACH(alignment::BinMetadata bin, bins)
    {
        if (!bin.isUnalignedBin())
        {
            bin.setIndex(ret.size());
            ret.push_back(bin);
        }
    }

    // put all the unaligned in the back if required
    if (keepUnaligned && putUnalignedInTheBack)
    {
        BOOST_FOREACH(alignment::BinMetadata bin, bins)
        {
            if (bin.isUnalignedBin())
            {
                breakUpBin(bin, partsCount, ret);
            }
        }
    }

    return ret;
}

Build::Build(const std::vector<std::string> &argv,
             const flowcell::FlowcellLayoutList &flowcellLayoutList,
             const flowcell::TileMetadataList &tileMetadataList,
             const flowcell::BarcodeMetadataList &barcodeMetadataList,
             const alignment::BinMetadataList &bins,
             const reference::SortedReferenceXmlList &sortedReferenceXmlList,
             const boost::filesystem::path outputDirectory,
             const unsigned maxLoaders,
             const unsigned maxComputers,
             const unsigned maxSavers,
             const build::GapRealignerMode realignGaps,
             const int bamGzipLevel,
             const double expectedBgzfCompressionRatio,
             const bool keepDuplicates,
             const bool clipSemialigned,
             const std::string &binRegexString,
             const bool keepUnknownAlignmentScore,
             const bool keepUnaligned,
             const bool putUnalignedInTheBack)
    :argv_(argv),
     tileMetadataList_(tileMetadataList),
     barcodeMetadataList_(barcodeMetadataList),
     bins_(breakUpUnalignedBin(filterBins(bins, binRegexString), maxComputers, keepUnaligned, putUnalignedInTheBack)),
     sortedReferenceXmlList_(sortedReferenceXmlList),
     outputDirectory_(outputDirectory),
     maxLoaders_(maxLoaders),
     maxComputers_(maxComputers),
     maxSavers_(maxSavers),
     bamGzipLevel_(bamGzipLevel),
     keepUnknownAlignmentScore_(keepUnknownAlignmentScore),
     keepDuplicates_(keepDuplicates),
     clipSemialigned_(clipSemialigned),
     realignGaps_(realignGaps),
     // assuming the last entry in the list contains the longest paths
     maxBinPathLength_(bins_.empty() ? 0 : std::max(bins_.back().getPathString().size(),
                                                    std::max(bins_.back().getFIdxFilePath().string().size(),
                                                             bins_.back().getRIdxFilePath().string().size()))),
     expectedBgzfCompressionRatio_(expectedBgzfCompressionRatio),
     threads_(maxComputers_ + maxLoaders_ + maxSavers_),
     contigList_(reference::loadContigs(sortedReferenceXmlList, threads_)),
     barcodeBamMapping_(mapBarcodesToFiles(outputDirectory_, barcodeMetadataList_)),
     bamFileStreams_(createOutputFileStreams(tileMetadataList_, barcodeMetadataList_)),
     // Stats must be initialized on the whole set of bins, so that bin index of each bin means something for it
     stats_(bins_, barcodeMetadataList_),
     maxReadLength_(getMaxReadLength(flowcellLayoutList)),
     threadBinSorters_(threads_.size()),
     threadBgzfBuffers_(threads_.size(), std::vector<std::vector<char> >(bamFileStreams_.size())),
     threadBgzfStreams_(threads_.size()),
     threadBamIndexParts_(threads_.size())
{
    computeSlotWaitingBins_.reserve(threads_.size());
    while(threadBgzfStreams_.size() < threads_.size())
    {
        threadBgzfStreams_.push_back(new boost::ptr_vector<boost::iostreams::filtering_ostream>(bamFileStreams_.size()));
    }
    while(threadBamIndexParts_.size() < threads_.size())
    {
        threadBamIndexParts_.push_back(new boost::ptr_vector<bam::BamIndexPart>(bamFileStreams_.size()));
    }
    threads_.execute(boost::bind(&Build::allocateThreadData, this, _1));

    // test if all the bins will fit in memory well
    {
        std::vector<std::vector<char> > testBgzfBuffers(bamFileStreams_.size());
        // no need to really block anything at the moment. just supply fake block to reserveBuffers so it compiles
        common::ScoopedMallocBlock fakeMallocBlock(common::ScoopedMallocBlock::Off);
        for(alignment::BinMetadataList::const_iterator binIterator = bins_.begin(); bins_.end() != binIterator; ++binIterator)
        {
            // NOTE: bin 0 contains unaligned data which can be large. As the unaligned reads are streamed directly from the source file
            // into compressed output, no need to check memory capacity for the bin 0.
            if (binIterator->getIndex())
            {
                boost::unique_lock<boost::mutex> lock(stateMutex_);
                const unsigned long bytesFailedToAllocate = reserveBuffers(lock, binIterator, bins_.end(), testBgzfBuffers, fakeMallocBlock, 0);
                if (bytesFailedToAllocate)
                {
                    BOOST_THROW_EXCEPTION(
                        common::MemoryException((
                            boost::format("%s requires %d bytes of memory for BAM generation.") %
                                *binIterator % bytesFailedToAllocate).str()));
                }
                BOOST_FOREACH(std::vector<char> &bgzfBuffer, testBgzfBuffers)
                {
                    // release memory as next bin test might fail if we don't
                    std::vector<char>().swap(bgzfBuffer);
                }
                threadBinSorters_.at(0).reset();
                boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams = threadBgzfStreams_.at(0);
                bgzfStreams.clear();
                boost::ptr_vector<bam::BamIndexPart> &bamIndexParts = threadBamIndexParts_.at(0);
                bamIndexParts.clear();

            }
        }
    }
}

void Build::allocateThreadData(const size_t threadNumber)
{
#ifdef HAVE_NUMA
    int runOnNode = threadNumber % (numa_max_node() + 1);
    ISAAC_ASSERT_MSG(8 * sizeof(unsigned long) >= unsigned(runOnNode), "numa node is too high: " << runOnNode);
    ISAAC_ASSERT_MSG(-1 != numa_run_on_node(runOnNode), "numa_run_on_node " << runOnNode <<
        " failed, errno: " << errno  << ":" << strerror(errno));
    unsigned long nodemask = 1UL << runOnNode;
    ISAAC_ASSERT_MSG(-1 != set_mempolicy(MPOL_BIND/*|MPOL_F_STATIC_NODES*/, &nodemask, sizeof(nodemask) * 8),
                     "set_mempolicy for nodemask: " << nodemask <<
                     " failed, errno: " << errno << ":" << strerror(errno));
#endif //HAVE_NUMA
}

void Build::run(common::ScoopedMallocBlock &mallocBlock)
{
    alignment::BinMetadataList::const_iterator nextUnprocessedBinIt(bins_.begin());
    alignment::BinMetadataList::const_iterator nextUnallocatedBinIt(bins_.begin());
    alignment::BinMetadataList::const_iterator nextUnloadedBinIt(bins_.begin());
    alignment::BinMetadataList::const_iterator nextUncompressedBinIt(bins_.begin());
    alignment::BinMetadataList::const_iterator nextUnsavedBinIt(bins_.begin());

    threads_.execute(boost::bind(&Build::sortBinParallel, this,
                                boost::ref(nextUnprocessedBinIt),
                                boost::ref(nextUnallocatedBinIt),
                                boost::ref(nextUnloadedBinIt),
                                boost::ref(nextUncompressedBinIt),
                                boost::ref(nextUnsavedBinIt),
                                bins_.end(),
                                boost::ref(mallocBlock),
                                _1));

    unsigned fileIndex = 0;
    BOOST_FOREACH(const boost::filesystem::path &bamFilePath, barcodeBamMapping_.getPaths())
    {
        // some of the streams are null_sink (that's when reference is unmapped for the sample).
        // this is the simplest way to ignore them...
        if (boost::filesystem::exists(bamFilePath))
        {
            boost::iostreams::filtering_ostream *stm = bamFileStreams_.at(fileIndex).get();
            if (stm)
            {
                bam::serializeBgzfFooter(*stm);
                ISAAC_THREAD_CERR << "BAM file generated: " << bamFilePath << "\n";
            }
        }
        ++fileIndex;
    }
}

void Build::dumpStats(const boost::filesystem::path &statsXmlPath)
{
    BuildStatsXml statsXml(sortedReferenceXmlList_, bins_, barcodeMetadataList_, stats_);
    std::ofstream os(statsXmlPath.string().c_str());
    if (!os) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: Unable to open file for writing: " + statsXmlPath.string()));
    }
    if (!(os << statsXml)) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: failed to store MatchFinder statistics in : " + statsXmlPath.string()));
    }
}

unsigned long Build::estimateOptimumFragmentsPerBin(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const unsigned long availableMemory,
    const double expectedBgzfCompressionRatio,
    const unsigned computeThreads)
{
//    const size_t maxFragmentIndexBytes = std::max(sizeof(io::RStrandOrShadowFragmentIndex),
//                                                         sizeof(io::FStrandFragmentIndex));

    const size_t maxFragmentBytes =
        io::FragmentHeader::getMinTotalLength(flowcell::getMaxReadLength(flowcellLayoutList));

    const size_t maxFragmentDedupedIndexBytes = sizeof(PackedFragmentBuffer::Index);
    const size_t maxFragmentCompressedBytes = maxFragmentBytes * expectedBgzfCompressionRatio;

    const size_t fragmentMemoryRequirements =
        // assume the initial indexes don't stay in memory for too long //maxFragmentIndexBytes              //index containing duplicates
        + maxFragmentBytes                 //data
        + maxFragmentDedupedIndexBytes     //deduplicated index
        + maxFragmentCompressedBytes       //bgzf chunk
        ;

    // reasonable amount of bins-in-progress to allow for no-delay input/compute/output overlap
//    const unsigned minOverlap = 3;;
    // try to increase granularity so that the CPU gets efficiently utilized.
    const unsigned minOverlap = computeThreads;
    return availableMemory / fragmentMemoryRequirements / minOverlap;
}

/**
 * \brief Attempts to reserve memory buffers required to process a bin.
 *
 * \return Non-zero size in bytes if the reservation failed.
 */
unsigned long Build::reserveBuffers(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataList::const_iterator thisThreadBinIt,
    const alignment::BinMetadataList::const_iterator binsEnd,
//    BinSorter &binSorter,
    std::vector<std::vector<char> > &bgzfBuffers,
    common::ScoopedMallocBlock &mallocBlock,
    const size_t threadNumber)
{
    common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
    boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams = threadBgzfStreams_.at(threadNumber);
    boost::ptr_vector<bam::BamIndexPart> &bamIndexParts = threadBamIndexParts_.at(threadNumber);
    const alignment::BinMetadata &bin = *thisThreadBinIt;
    common::ScoopedMallocBlockUnblock unblockMalloc(mallocBlock);
    try
    {
        threadBinSorters_.at(threadNumber) = boost::shared_ptr<BinSorter>(
            new BinSorter(keepDuplicates_, clipSemialigned_, barcodeBamMapping_, tileMetadataList_, barcodeMetadataList_,
                          maxReadLength_, realignGaps_, contigList_, keepUnknownAlignmentScore_,
                          bin));

        threadBinSorters_.at(threadNumber)->reservePathBuffers(maxBinPathLength_);

        BOOST_FOREACH(std::vector<char> &bgzfBuffer, bgzfBuffers)
        {
            const unsigned outputFileIndex = &bgzfBuffer - &bgzfBuffers.front();
            bgzfBuffer.reserve(estimateBinCompressedDataRequirements(bin, outputFileIndex));
        }

        ISAAC_ASSERT_MSG(!bgzfStreams.size(), "Expecting empty pool of streams");
        while(bgzfStreams.size() < bamFileStreams_.size())
        {
            bgzfStreams.push_back(new boost::iostreams::filtering_ostream);
            bgzfStreams.back().push(bam::BgzfCompressor(bamGzipLevel_));
            bgzfStreams.back().push(
                boost::iostreams::back_insert_device<std::vector<char> >(
                    threadBgzfBuffers_.at(threadNumber).at(bgzfStreams.size()-1)));
//            bgzfStreams.back().push(boost::iostreams::null_sink());
        }

        ISAAC_ASSERT_MSG(!bamIndexParts.size(), "Expecting empty pool of bam index parts");
        while(bamIndexParts.size() < bamFileStreams_.size())
        {
            bamIndexParts.push_back(new bam::BamIndexPart);
        }



    }
    catch(std::bad_alloc &e)
    {
        errno = 0;
        bgzfStreams.clear();
        bamIndexParts.clear();
        // give a chance other threads to allocate what they need... TODO: this is not required anymore as allocation happens orderly
        threadBinSorters_.at(threadNumber).reset();
        unsigned long totalBuffersNeeded = 0UL;
        BOOST_FOREACH(std::vector<char> &bgzfBuffer, bgzfBuffers)
        {
            std::vector<char>().swap(bgzfBuffer);
            const unsigned outputFileIndex = &bgzfBuffer - &bgzfBuffers.front();
            totalBuffersNeeded += estimateBinCompressedDataRequirements(bin, outputFileIndex);
        }
        return BinSorter::getMemoryRequirements(bin) + totalBuffersNeeded;
    }
    return 0;
}


void Build::allocateBin(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataList::const_iterator thisThreadBinIt,
    const alignment::BinMetadataList::const_iterator binsEnd,
    alignment::BinMetadataList::const_iterator &nextUnallocatedBinIt,
    std::vector<std::vector<char> > &bgzfBuffers,
    common::ScoopedMallocBlock &mallocBlock,
    const size_t threadNumber)
{
    unsigned long requiredMemory = 0;
    while(nextUnallocatedBinIt != thisThreadBinIt ||
          0 != (requiredMemory = reserveBuffers(lock, thisThreadBinIt, binsEnd, bgzfBuffers, mallocBlock, threadNumber)))
    {
        if (nextUnallocatedBinIt == thisThreadBinIt)
        {
            ISAAC_THREAD_CERR << "WARNING: Holding up processing of bin: " <<
                thisThreadBinIt->getPath() << " until " << requiredMemory <<
                " bytes of allowed memory is available." << std::endl;
        }
        stateChangedCondition_.wait(lock);
    }
    ++nextUnallocatedBinIt;
    stateChangedCondition_.notify_all();

}

void Build::waitForLoadSlot(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataList::const_iterator thisThreadBinIt,
    const alignment::BinMetadataList::const_iterator binsEnd,
    alignment::BinMetadataList::const_iterator &nextUnloadedBinIt,
    std::vector<std::vector<char> > &bgzfBuffers,
    common::ScoopedMallocBlock &mallocBlock,
    const size_t threadNumber)
{
    while(nextUnloadedBinIt != thisThreadBinIt || !maxLoaders_)
    {
        if (nextUnloadedBinIt == thisThreadBinIt)
        {
            ISAAC_THREAD_CERR << "WARNING: Holding up processing of bin: " <<
                thisThreadBinIt->getPath() << " until a load slot is available" << std::endl;
        }
        stateChangedCondition_.wait(lock);
    }

    ++nextUnloadedBinIt;
    --maxLoaders_;
    stateChangedCondition_.notify_all();
}

void Build::returnLoadSlot()
{
    ++maxLoaders_;
    stateChangedCondition_.notify_all();
}

void Build::waitForComputeSlot(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataList::const_iterator thisThreadBinIt,
    alignment::BinMetadataList::const_iterator &nextUncompressedBinIt)
{
    const unsigned binIndex = thisThreadBinIt->getIndex();
    computeSlotWaitingBins_.push_back(binIndex);

    while(!maxComputers_ ||
        binIndex != *std::min_element(computeSlotWaitingBins_.begin(), computeSlotWaitingBins_.end()))
    {
        if (nextUncompressedBinIt == thisThreadBinIt)
        {
            ISAAC_THREAD_CERR << "WARNING: Holding up processing of bin: " <<
                thisThreadBinIt->getPath() << " until a compute slot is available." << std::endl;
        }

        stateChangedCondition_.wait(lock);
    }
    --maxComputers_;
    computeSlotWaitingBins_.erase(std::find(computeSlotWaitingBins_.begin() ,computeSlotWaitingBins_.end(), binIndex));
    ++nextUncompressedBinIt;
    stateChangedCondition_.notify_all();
}

void Build::returnComputeSlot()
{
    ++maxComputers_;
    stateChangedCondition_.notify_all();
}

void Build::waitForSaveSlot(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataList::const_iterator thisThreadBinIt,
    alignment::BinMetadataList::const_iterator &nextUnsavedBinIt)
{
    while(nextUnsavedBinIt != thisThreadBinIt)
    {
        stateChangedCondition_.wait(lock);
    }
}

void Build::returnSaveSlot(
    alignment::BinMetadataList::const_iterator &nextUnsavedBinIt)
{
    ++nextUnsavedBinIt;
    stateChangedCondition_.notify_all();
}

void Build::sortBinParallel(alignment::BinMetadataList::const_iterator &nextUnprocessedBinIt,
                            alignment::BinMetadataList::const_iterator &nextUnallocatedBinIt,
                            alignment::BinMetadataList::const_iterator &nextUnloadedBinIt,
                            alignment::BinMetadataList::const_iterator &nextUncompressedBinIt,
                            alignment::BinMetadataList::const_iterator &nextUnsavedBinIt,
                            const alignment::BinMetadataList::const_iterator binsEnd,
                            common::ScoopedMallocBlock &mallocBlock,
                            const size_t threadNumber)
{
    boost::unique_lock<boost::mutex> lock(stateMutex_);
    while(binsEnd != nextUnprocessedBinIt)
    {
        alignment::BinMetadataList::const_iterator thisThreadBinIt;
        thisThreadBinIt = nextUnprocessedBinIt++;

        std::vector<std::vector<char> > &bgzfBuffers = threadBgzfBuffers_.at(threadNumber);

        // wait and allocate memory required for loading and compressing this bin
        allocateBin(lock, thisThreadBinIt, binsEnd, nextUnallocatedBinIt, bgzfBuffers, mallocBlock, threadNumber);
        BinSorter &indexedBin(*threadBinSorters_.at(threadNumber));
        waitForLoadSlot(lock, thisThreadBinIt, binsEnd, nextUnloadedBinIt, bgzfBuffers, mallocBlock, threadNumber);
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&Build::returnLoadSlot, this))
        {
            {
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                indexedBin.load();
            }
        }

        waitForComputeSlot(lock, thisThreadBinIt, nextUncompressedBinIt);
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&Build::returnComputeSlot, this))
        {
            {
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                processBin(indexedBin, threadNumber);
                boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams = threadBgzfStreams_.at(threadNumber);
                bgzfStreams.clear();
            }
            // give back some memory to allow other threads to load
            // data while we're waiting for our turn to save
            threadBinSorters_.at(threadNumber).reset();
        }

        // wait for our turn to store bam data
        waitForSaveSlot(lock, thisThreadBinIt, nextUnsavedBinIt);
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&Build::returnSaveSlot, this, boost::ref(nextUnsavedBinIt)))
        {
            saveAndReleaseBuffers(lock, bgzfBuffers, threadBamIndexParts_.at(threadNumber), thisThreadBinIt->getPath());
        }
    }
}

unsigned long Build::processBin(
    BinSorter &indexedBin,
    const unsigned threadNumber)
{
    const unsigned long unique = indexedBin.process(stats_);
    if (unique)
    {
        indexedBin.reorderForBam();
        indexedBin.serialize(threadBgzfStreams_.at(threadNumber),
                             threadBamIndexParts_.at(threadNumber));
    }
    return unique;
}

/**
 * \brief Save bgzf compressed buffers into corresponding sample files and and release associated memory
 */
void Build::saveAndReleaseBuffers(
    boost::unique_lock<boost::mutex> &lock,
    std::vector<std::vector<char> > &bgzfBuffers,
    boost::ptr_vector<bam::BamIndexPart> &bamIndexParts,
    const boost::filesystem::path &filePath)
{
    BOOST_FOREACH(std::vector<char> &bgzfBuffer, bgzfBuffers)
    {
        {
            common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
            unsigned index = &bgzfBuffer - &bgzfBuffers.front();
            boost::iostreams::filtering_ostream *stm = bamFileStreams_.at(index).get();
            if (!stm)
            {
                ISAAC_ASSERT_MSG(bgzfBuffer.empty(), "Unexpected data for bam file belonging to a sample with unmapped reference");
            }
            else
            {
                saveBuffer(bgzfBuffer, *stm, bamIndexParts.at(index), bamIndexes_.at(index), filePath);
            }
        }
        // release rest of the memory that was reserved for this bin
        std::vector<char>().swap(bgzfBuffer);
    }
    bamIndexParts.clear();
}

void Build::saveBuffer(
    const std::vector<char> &bgzfBuffer,
    boost::iostreams::filtering_ostream &bamStream,
    const bam::BamIndexPart &bamIndexPart,
    bam::BamIndex &bamIndex,
    const boost::filesystem::path &filePath)
{
    ISAAC_THREAD_CERR << "Saving " << bgzfBuffer.size() << " bytes of sorted data for bin " << filePath << std::endl;
    const clock_t start = clock();
    if(!bamStream.write(&bgzfBuffer.front(), bgzfBuffer.size())/* ||
        !bamStream.strict_sync()*/){
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Failed to write bgzf block of %d bytes into bam stream") % bgzfBuffer.size()).str()));
    }
    bamIndex.processIndexPart( bamIndexPart, bgzfBuffer );

    ISAAC_THREAD_CERR << "Saving " << bgzfBuffer.size() << " bytes of sorted data for bin " << filePath << " done in " << (clock() - start) / 1000 << "ms\n";
}

} // namespace build
} // namespace isaac
