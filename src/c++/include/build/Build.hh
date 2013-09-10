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
#include "alignment/TemplateLengthStatistics.hh"
#include "build/BarcodeBamMapping.hh"
#include "build/BinSorter.hh"
#include "build/BuildStats.hh"
#include "build/BuildContigMap.hh"
#include "common/Threads.hpp"
#include "reference/SortedReferenceMetadata.hh"


namespace isaac
{
namespace build
{

class Build
{
    const std::vector<std::string> &argv_;
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    const flowcell::TileMetadataList &tileMetadataList_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    alignment::BinMetadataList unalignedBinParts_;
    const alignment::BinMetadataCRefList bins_;
    const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics_;
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList_;
    const BuildContigMap contigMap_;
    const boost::filesystem::path outputDirectory_;
    unsigned maxLoaders_;
    unsigned maxComputers_;
    std::vector<unsigned> computeSlotWaitingBins_;
    const unsigned maxSavers_;
    const int bamGzipLevel_;
    const std::vector<std::string> &bamHeaderTags_;
    // forcedDodgyAlignmentScore_ gets assigned to reads that have their scores at ushort -1
    const unsigned char forcedDodgyAlignmentScore_;
    const bool singleLibrarySamples_;
    const bool keepDuplicates_;
    const bool markDuplicates_;
    const bool realignGapsVigorously_;
    const bool realignDodgyFragments_;
    const unsigned realignedGapsPerFragment_;
    const bool clipSemialigned_;
    const build::GapRealignerMode realignGaps_;
    const double expectedBgzfCompressionRatio_;
    const unsigned maxReadLength_;
    const IncludeTags includeTags_;
    const bool pessimisticMapQ_;

    boost::mutex stateMutex_;
    boost::condition_variable stateChangedCondition_;

    common::ThreadVector threads_;

    const std::vector<std::vector<reference::Contig> > contigList_;
    //pair<[barcode], [output file]>, first maps barcode indexes to unique paths in second
    BarcodeBamMapping barcodeBamMapping_;
    //[output file], one stream per bam file path
    boost::ptr_vector<bam::BamIndex> bamIndexes_;
    std::vector<boost::shared_ptr<std::ofstream> > bamFileStreams_;

    BuildStats stats_;

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
          const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
          const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
          const boost::filesystem::path outputDirectory,
          const unsigned maxLoaders,
          const unsigned maxComputers,
          const unsigned maxSavers,
          const build::GapRealignerMode realignGaps,
          const int bamGzipLevel,
          const std::vector<std::string> &bamHeaderTags,
          const double expectedBgzfCompressionRatio,
          const bool singleLibrarySamples,
          const bool keepDuplicates,
          const bool markDuplicates,
          const bool realignGapsVigorously,
          const bool realignDodgyFragments,
          const unsigned realignedGapsPerFragment,
          const bool clipSemialigned,
          const std::string &binRegexString,
          const unsigned char forcedDodgyAlignmentScore,
          const bool keepUnaligned,
          const bool putUnalignedInTheBack,
          const IncludeTags includeTags,
          const bool pessimisticMapQ);

    void run(common::ScoopedMallocBlock &mallocBlock);

    void dumpStats(const boost::filesystem::path &statsXmlPath);

    static unsigned long estimateOptimumFragmentsPerBin(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned long availableMemory,
        const double expectedBgzfCompressionRatio,
        const unsigned computeThreads);

    const BarcodeBamMapping &getBarcodeBamMapping() const {return barcodeBamMapping_;}
private:
    std::vector<boost::shared_ptr<std::ofstream> >  createOutputFileStreams(
        const flowcell::TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList);

    unsigned long reserveBuffers(
        boost::unique_lock<boost::mutex> &lock,
        const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
        const alignment::BinMetadataCRefList::const_iterator binsEnd,
        common::ScoopedMallocBlock &mallocBlock,
        const size_t threadNumber);

    void allocateThreadData(const size_t threadNumber);

    void allocateBin(
        boost::unique_lock<boost::mutex> &lock,
        const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
        const alignment::BinMetadataCRefList::const_iterator binsEnd,
        alignment::BinMetadataCRefList::const_iterator &nextUnallocatedBinIt,
        common::ScoopedMallocBlock &mallocBlock,
        const size_t threadNumber);

    void waitForLoadSlot(
        boost::unique_lock<boost::mutex> &lock,
        const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
        const alignment::BinMetadataCRefList::const_iterator binsEnd,
        alignment::BinMetadataCRefList::const_iterator &nextUnloadedBinIt,
        common::ScoopedMallocBlock &mallocBlock,
        const size_t threadNumber);

    void returnLoadSlot();

    void waitForComputeSlot(
        boost::unique_lock<boost::mutex> &lock,
        const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
        alignment::BinMetadataCRefList::const_iterator &nextUncompressedBinIt);

    void returnComputeSlot();

    void waitForSaveSlot(
        boost::unique_lock<boost::mutex> &lock,
        const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
        alignment::BinMetadataCRefList::const_iterator &nextUnserializedBinIt);

    void returnSaveSlot(
        alignment::BinMetadataCRefList::const_iterator &nextUnserializedBinIt);

    void sortBinParallel(alignment::BinMetadataCRefList::const_iterator &nextUnprocessedBinIt,
                         alignment::BinMetadataCRefList::const_iterator &nextUnallocatedBinIt,
                         alignment::BinMetadataCRefList::const_iterator &nextUnloadedBinIt,
                         alignment::BinMetadataCRefList::const_iterator &nextUncompressedBinIt,
                         alignment::BinMetadataCRefList::const_iterator &nextUnserializedBinIt,
                         common::ScoopedMallocBlock &mallocBlock,
                         const size_t threadNumber);

    unsigned long processBin(
        BinSorter &indexedBin,
        const unsigned threadNumber);

    void saveAndReleaseBuffers(
        boost::unique_lock<boost::mutex> &lock,
        const boost::filesystem::path &filePath,
        const size_t threadNumber);

    void saveBuffer(
        const std::vector<char> &bgzfBuffer,
        std::ostream &bamStream,
        const bam::BamIndexPart &bamIndexPart,
        bam::BamIndex &bamIndex,
        const boost::filesystem::path &filePath);

    unsigned long estimateBinCompressedDataRequirements(
        const alignment::BinMetadata & binMetadata,
        const unsigned outputFileIndex) const;

    void testBinsFitInRam();
};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BUILD_HH
