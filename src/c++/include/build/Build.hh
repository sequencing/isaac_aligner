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
 ** \file Build.hh
 **
 ** Reorders alingments and stores them in bam file.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BUILD_HH
#define iSAAC_BUILD_BUILD_HH

#include <boost/filesystem.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/thread.hpp>

#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/TileMetadata.hh"
#include "alignment/BinMetadata.hh"
#include "build/BarcodeBamMapping.hh"
#include "build/BinSorter.hh"
#include "build/BuildStats.hh"
#include "common/Threads.hpp"
#include "reference/SortedReferenceXml.hh"


namespace isaac
{
namespace build
{

namespace bfs = boost::filesystem;

class Build
{
    const std::vector<std::string> &argv_;
    const flowcell::TileMetadataList &tileMetadataList_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const alignment::BinMetadataList bins_;
    const reference::SortedReferenceXmlList &sortedReferenceXmlList_;
    const boost::filesystem::path outputDirectory_;
    unsigned maxLoaders_;
    unsigned maxComputers_;
    std::vector<unsigned> computeSlotWaitingBins_;
    const unsigned maxSavers_;
    const int bamGzipLevel_;
    const bool keepUnknownAlignmentScore_;
    const bool keepDuplicates_;
    const bool clipSemialigned_;
    const build::GapRealignerMode realignGaps_;
    const size_t maxBinPathLength_;
    const double expectedBgzfCompressionRatio_;

    boost::mutex stateMutex_;
    boost::condition_variable stateChangedCondition_;

    common::ThreadVector threads_;

    const std::vector<std::vector<reference::Contig> > contigList_;
    //pair<[barcode], [output file]>, first maps barcode indexes to unique paths in second
    BarcodeBamMapping barcodeBamMapping_;
    //[output file], one stream per bam file path
    boost::ptr_vector<bam::BamIndex> bamIndexes_;
    std::vector<boost::shared_ptr<boost::iostreams::filtering_ostream> > bamFileStreams_;

    BuildStats stats_;

    const unsigned maxReadLength_;
    //[thread]
    std::vector<boost::shared_ptr<BinSorter> > threadBinSorters_;
    //[thread][bam file][byte]
    std::vector<std::vector<std::vector<char> > > threadBgzfBuffers_;
    // Geometry: [thread][bam file]. Streams for compressing bam data into threadBgzfBuffers_
    boost::ptr_vector<boost::ptr_vector<boost::iostreams::filtering_ostream> > threadBgzfStreams_;
    boost::ptr_vector<boost::ptr_vector<bam::BamIndexPart> > threadBamIndexParts_;

public:
    Build(const std::vector<std::string> &argv,
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
          const bool putUnalignedInTheBack);

    void run(common::ScoopedMallocBlock &mallocBlock);

    void dumpStats(const boost::filesystem::path &statsXmlPath);

    static unsigned long estimateOptimumFragmentsPerBin(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned long availableMemory,
        const double expectedBgzfCompressionRatio,
        const unsigned computeThreads);

    const BarcodeBamMapping &getBarcodeBamMapping() const {return barcodeBamMapping_;}
private:
    std::vector<boost::shared_ptr<boost::iostreams::filtering_ostream> >  createOutputFileStreams(
        const flowcell::TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList);

    unsigned long reserveBuffers(
        boost::unique_lock<boost::mutex> &lock,
        const alignment::BinMetadataList::const_iterator thisThreadBinIt,
        const alignment::BinMetadataList::const_iterator binsEnd,
//        BinSorter &binSorter,
        std::vector<std::vector<char> > &bgzfBuffers,
        common::ScoopedMallocBlock &mallocBlock,
        const size_t threadNumber);

    void allocateThreadData(const size_t threadNumber);

    void allocateBin(
        boost::unique_lock<boost::mutex> &lock,
        const alignment::BinMetadataList::const_iterator thisThreadBinIt,
        const alignment::BinMetadataList::const_iterator binsEnd,
        alignment::BinMetadataList::const_iterator &nextUnallocatedBinIt,
        std::vector<std::vector<char> > &bgzfBuffers,
        common::ScoopedMallocBlock &mallocBlock,
        const size_t threadNumber);

    void waitForLoadSlot(
        boost::unique_lock<boost::mutex> &lock,
        const alignment::BinMetadataList::const_iterator thisThreadBinIt,
        const alignment::BinMetadataList::const_iterator binsEnd,
        alignment::BinMetadataList::const_iterator &nextUnloadedBinIt,
        std::vector<std::vector<char> > &bgzfBuffers,
        common::ScoopedMallocBlock &mallocBlock,
        const size_t threadNumber);

    void returnLoadSlot();

    void waitForComputeSlot(
        boost::unique_lock<boost::mutex> &lock,
        const alignment::BinMetadataList::const_iterator thisThreadBinIt,
        alignment::BinMetadataList::const_iterator &nextUncompressedBinIt);

    void returnComputeSlot();

    void waitForSaveSlot(
        boost::unique_lock<boost::mutex> &lock,
        const alignment::BinMetadataList::const_iterator thisThreadBinIt,
        alignment::BinMetadataList::const_iterator &nextUnserializedBinIt);

    void returnSaveSlot(
        alignment::BinMetadataList::const_iterator &nextUnserializedBinIt);

    void sortBinParallel(alignment::BinMetadataList::const_iterator &nextUnprocessedBinIt,
                         alignment::BinMetadataList::const_iterator &nextUnallocatedBinIt,
                         alignment::BinMetadataList::const_iterator &nextUnloadedBinIt,
                         alignment::BinMetadataList::const_iterator &nextUncompressedBinIt,
                         alignment::BinMetadataList::const_iterator &nextUnserializedBinIt,
                         const alignment::BinMetadataList::const_iterator binsEnd,
                         common::ScoopedMallocBlock &mallocBlock,
                         const size_t threadNumber);

    unsigned long processBin(
        BinSorter &indexedBin,
        const unsigned threadNumber);

    void saveAndReleaseBuffers(
        boost::unique_lock<boost::mutex> &lock,
        std::vector<std::vector<char> > &bgzfBuffers,
        boost::ptr_vector<bam::BamIndexPart> &bamIndexParts,
        const boost::filesystem::path &filePath);

    void saveBuffer(
        const std::vector<char> &bgzfBuffer,
        boost::iostreams::filtering_ostream &bamStream,
        const bam::BamIndexPart &bamIndexPart,
        bam::BamIndex &bamIndex,
        const boost::filesystem::path &filePath);

    unsigned long estimateBinCompressedDataRequirements(
        const alignment::BinMetadata & binMetadata,
        const unsigned outputFileIndex) const;
};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BUILD_HH
