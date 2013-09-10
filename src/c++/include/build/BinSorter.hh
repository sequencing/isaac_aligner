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
 ** \file BinSorter.hh
 **
 ** Performs sorting and duplicate marking on a single alignment bin.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BIN_SORTER_HH
#define iSAAC_BUILD_BIN_SORTER_HH

#include <numeric>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "alignment/BinMetadata.hh"
#include "build/BamSerializer.hh"
#include "build/DuplicatePairEndFilter.hh"
#include "build/FragmentIndex.hh"
#include "build/GapRealigner.hh"
#include "build/NotAFilter.hh"
#include "build/DuplicateFragmentIndexFiltering.hh"
#include "build/PackedFragmentBuffer.hh"
#include "io/FileBufCache.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace build
{

enum GapRealignerMode
{
    /// don't realign
    REALIGN_NONE,
    /// Realign against gaps found within the sample
    REALIGN_SAMPLE,
    /// Realign against gaps found in all samples of the same project
    REALIGN_PROJECT,
    /// Realign against all gaps present in the data
    REALIGN_ALL
};

class BinSorter : std::vector<PackedFragmentBuffer::Index>
{
    typedef std::vector<PackedFragmentBuffer::Index> BaseType;

public:
    BinSorter(
        const bool singleLibrarySamples,
        const bool keepDuplicates,
        const bool markDuplicates,
        const bool realignGapsVigorously,
        const bool realignDodgyFragments,
        const unsigned realignedGapsPerFragment,
        const bool clipSemialigned,
        const BarcodeBamMapping &barcodeBamMapping,
        const flowcell::TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        const BuildContigMap &contigMap,
        const unsigned maxReadLength,
        const build::GapRealignerMode realignGaps,
        const std::vector<std::vector<isaac::reference::Contig> > &contigList,
        const unsigned char forcedDodgyAlignmentScore,
        const alignment::BinMetadata &bin,
        const unsigned binStatsIndex,
        const flowcell::FlowcellLayoutList &flowCellLayoutList,
        const IncludeTags includeTags,
        const bool pessimisticMapQ) :
            singleLibrarySamples_(singleLibrarySamples),
            keepDuplicates_(keepDuplicates),
            markDuplicates_(markDuplicates),
            bin_(bin),
            binStatsIndex_(binStatsIndex),
            barcodeBamMapping_(barcodeBamMapping),
            bamSerializer_(barcodeBamMapping_.getSampleIndexMap(), tileMetadataList, barcodeMetadataList,
                           contigMap,
                           maxReadLength, forcedDodgyAlignmentScore, flowCellLayoutList, includeTags, pessimisticMapQ),
            fileBuf_(1, std::ios_base::binary|std::ios_base::in),
            realignGaps_(realignGaps),
            realignerGaps_(getGapGroupsCount()),
            gapRealigner_(
                realignGapsVigorously, realignDodgyFragments, realignedGapsPerFragment, 3, 4, 0, clipSemialigned,
                barcodeMetadataList, barcodeTemplateLengthStatistics, contigList),
            dataDistribution_(bin_.getDataDistribution())
    {
        data_.resize(bin_);

        std::vector<PackedFragmentBuffer::Index>::reserve(bin_.getTotalElements());
        seIdxFileContent_.reserve(bin_.getSeIdxElements());
        rIdxFileContent_.reserve(bin_.getRIdxElements());
        fIdxFileContent_.reserve(bin_.getFIdxElements());
        if (REALIGN_NONE != realignGaps_)
        {
            gapRealigner_.reserve(bin_);
            reserveGaps(bin_, barcodeMetadataList);

        }
        fileBuf_.reservePathBuffers(bin_.getPathString().size());
    }

    void reserveGaps(
        const alignment::BinMetadata& bin,
        const flowcell::BarcodeMetadataList &barcodeMetadataList);

    static unsigned long getMemoryRequirements(const alignment::BinMetadata& bin)
    {
        return PackedFragmentBuffer::getMemoryRequirements(bin) +
            bin.getSeIdxElements() * sizeof(SeFragmentIndex) +
            bin.getRIdxElements() * sizeof(RStrandOrShadowFragmentIndex) +
            bin.getFIdxElements() * sizeof(FStrandFragmentIndex) +
            bin.getTotalElements() * sizeof(PackedFragmentBuffer::Index);
    }

    void unreserveIndexes()
    {
        std::vector<SeFragmentIndex>().swap(seIdxFileContent_);
        std::vector<RStrandOrShadowFragmentIndex>().swap(rIdxFileContent_);
        std::vector<FStrandFragmentIndex>().swap(fIdxFileContent_);
    }

    void load()
    {
        ISAAC_THREAD_CERR << "Loading unsorted data" << std::endl;
        const clock_t startLoad = clock();

        loadData();

        ISAAC_THREAD_CERR << "Loading unsorted data done in " << (clock() - startLoad) / 1000 << "ms" << std::endl;
    }

    void reorderForBam()
    {
        ISAAC_THREAD_CERR << "Sorting offsets" << std::endl;
        if (REALIGN_NONE != realignGaps_)
        {
            // update all index record positions before sorting as they may have been messed up by gap realignment
            BOOST_FOREACH(PackedFragmentBuffer::Index &index, std::make_pair(indexBegin(), indexEnd()))
            {
                io::FragmentAccessor &fragment = data_.getFragment(index);
                index.pos_ = fragment.fStrandPosition_;
            }
        }
        const clock_t startSortOffsets = clock();
        std::sort(begin(), end(), boost::bind(&PackedFragmentBuffer::orderForBam, boost::ref(data_), _1, _2));
        ISAAC_THREAD_CERR << "Sorting offsets" << " done in " << (clock() - startSortOffsets) / 1000 << "ms" << std::endl;
    }

    unsigned long process(
        BuildStats &buildStats)
    {
        resolveDuplicates(buildStats);
        unreserveIndexes();
        if (!isUnalignedBin() && REALIGN_NONE != realignGaps_)
        {
            collectGaps();
            realignGaps();
        }

        return getUniqueRecordsCount();
    }

    unsigned long serialize(boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams,
                            boost::ptr_vector<bam::BamIndexPart> &bamIndexParts);

    unsigned getBinIndex() const
    {
        return bin_.getIndex();
    }
private:
    const bool singleLibrarySamples_;
    const bool keepDuplicates_;
    const bool markDuplicates_;
    const alignment::BinMetadata &bin_;
    const unsigned binStatsIndex_;
    const BarcodeBamMapping &barcodeBamMapping_;
    BamSerializer bamSerializer_;
    std::vector<SeFragmentIndex> seIdxFileContent_;
    std::vector<RStrandOrShadowFragmentIndex> rIdxFileContent_;
    std::vector<FStrandFragmentIndex> fIdxFileContent_;
    PackedFragmentBuffer data_;
    io::FileBufCache<io::FileBufWithReopen> fileBuf_;
    const GapRealignerMode realignGaps_;
    std::vector<RealignerGaps> realignerGaps_;
    GapRealigner gapRealigner_;
    alignment::BinDataDistribution dataDistribution_;

    void loadData();
    void loadUnalignedData();
    void loadAlignedData();
    const io::FragmentAccessor &loadFragment(std::istream &isData, unsigned long &offset);
    bool isUnalignedBin() const {return bin_.isUnalignedBin();}
    unsigned long getUniqueRecordsCount() const {return isUnalignedBin() ? bin_.getTotalElements() : size();}

    void resolveDuplicates(BuildStats &buildStats);
    void collectGaps();
    void realignGaps();

    BaseType::iterator indexBegin() {return begin();}
    BaseType::iterator indexEnd() {return end();}

    unsigned getGapGroupIndex(const unsigned barcode) const;
    unsigned getGapGroupsCount() const;

};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BIN_SORTER_HH
