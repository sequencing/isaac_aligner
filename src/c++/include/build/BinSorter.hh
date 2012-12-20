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
#include <boost/ptr_container/ptr_vector.hpp>

#include "alignment/BinMetadata.hh"
#include "build/BamSerializer.hh"
#include "build/DuplicatePairEndFilter.hh"
#include "build/GapRealigner.hh"
#include "build/SingleEndFilter.hh"
#include "build/DuplicateFragmentIndexFiltering.hh"
#include "build/PackedFragmentBuffer.hh"
#include "io/FileBufCache.hh"
#include "io/FragmentIndex.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace build
{

enum GapRealignerMode
{
    /// don't realign
    GAP_REALIGNER_OFF,
    /// brief Perform full realignment discarding all the existing gaps
    GAP_REALIGNER_ON
};

class BinSorter : std::vector<PackedFragmentBuffer::Index>
{
    typedef std::vector<PackedFragmentBuffer::Index> BaseType;

public:
    BinSorter(
        const bool keepDuplicates,
        const bool clipSemialigned,
        const BarcodeBamMapping &barcodeBamMapping,
        const flowcell::TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const unsigned maxReadLength,
        const build::GapRealignerMode realignGaps,
        const std::vector<std::vector<reference::Contig> > &contigList,
        const bool keepUnknownAlignmentScore,
        const alignment::BinMetadata &bin) :
            keepDuplicates_(keepDuplicates),
            bin_(bin),
            bamSerializer_(barcodeBamMapping.getIndexMap(), tileMetadataList, barcodeMetadataList,
                           maxReadLength, keepUnknownAlignmentScore),
            fileBuf_(1, std::ios_base::binary|std::ios_base::in),
            realignGaps_(realignGaps),
            gapRealigner_(clipSemialigned, barcodeMetadataList, contigList, barcodeBamMapping),
            dataDistribution_(barcodeMetadataList.size(), bin.getLength())
    {
        data_.resize(bin_);

        std::vector<PackedFragmentBuffer::Index>::reserve(bin_.getTotalElements());
        seIdxFileContent_.reserve(bin_.getSeIdxElements());
        rIdxFileContent_.reserve(bin_.getRIdxElements());
        fIdxFileContent_.reserve(bin_.getFIdxElements());
        if (GAP_REALIGNER_OFF != realignGaps_)
        {
            gapRealigner_.reserve(bin_);
        }
        dataDistribution_.reserve(bin_.getDataDistribution().size());
    }

    void reservePathBuffers(const size_t maxBinPathLength)
    {
        fileBuf_.reservePathBuffers(maxBinPathLength);
    }

    static unsigned long getMemoryRequirements(const alignment::BinMetadata& bin)
    {
        return PackedFragmentBuffer::getMemoryRequirements(bin) +
            bin.getSeIdxElements() * sizeof(io::SeFragmentIndex) +
            bin.getRIdxElements() * sizeof(io::RStrandOrShadowFragmentIndex) +
            bin.getFIdxElements() * sizeof(io::FStrandFragmentIndex) +
            bin.getTotalElements() * sizeof(PackedFragmentBuffer::Index);
    }

    void unreserve()
    {
        data_.unreserve();
        unreserveIndexes();
        std::vector<PackedFragmentBuffer::Index>().swap(*this);
        gapRealigner_.unreserve();
    }

    void unreserveIndexes()
    {
        std::vector<io::SeFragmentIndex>().swap(seIdxFileContent_);
        std::vector<io::RStrandOrShadowFragmentIndex>().swap(rIdxFileContent_);
        std::vector<io::FStrandFragmentIndex>().swap(fIdxFileContent_);
    }

    void load()
    {
        ISAAC_THREAD_CERR << "Loading unsorted data" << std::endl;
        const clock_t startLoad = clock();

        // indexes are loaded before data as they need their offsets patched during data loading
        loadIndex(bin_.getSeIdxFilePath(), bin_.getSeIdxElements(), seIdxFileContent_);
        loadIndex(bin_.getRIdxFilePath(), bin_.getRIdxElements(), rIdxFileContent_);
        loadIndex(bin_.getFIdxFilePath(), bin_.getFIdxElements(), fIdxFileContent_);
        loadData();

        ISAAC_THREAD_CERR << "Loading unsorted data done in " << (clock() - startLoad) / 1000 << "ms" << std::endl;
    }

    void reorderForBam()
    {
        ISAAC_THREAD_CERR << "Sorting offsets" << std::endl;
        if (GAP_REALIGNER_OFF != realignGaps_)
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
        SingleEndFilter().filterInput(data_, seIdxFileContent_.begin(), seIdxFileContent_.end(), buildStats, bin_.getIndex(), std::back_inserter<BaseType>(*this));
        DuplicatePairEndFilter(keepDuplicates_).filterInput(data_, rIdxFileContent_.begin(), rIdxFileContent_.end(), buildStats, bin_.getIndex(), std::back_inserter<BaseType>(*this));
        DuplicatePairEndFilter(keepDuplicates_).filterInput(data_, fIdxFileContent_.begin(), fIdxFileContent_.end(), buildStats, bin_.getIndex(), std::back_inserter<BaseType>(*this));

        unreserveIndexes();

        if (!isUnalignedBin() && GAP_REALIGNER_ON == realignGaps_)
        {
            ISAAC_THREAD_CERR << "Finalizing " << gapRealigner_.getTotalGapsCount() << " gaps."  << std::endl;
            for(std::vector<char>::const_iterator p = data_.begin(); p != data_.end();)
            {
                const io::FragmentAccessor &fragment = *reinterpret_cast<const io::FragmentAccessor *>(&*p);
                gapRealigner_.addGapsFromFragment(fragment);
                p += fragment.getTotalLength();
            }

            gapRealigner_.finalizeGaps();

            ISAAC_THREAD_CERR << "Realigning against " << gapRealigner_.getTotalGapsCount() << " unique gaps."  << std::endl;
            BOOST_FOREACH(PackedFragmentBuffer::Index &index, std::make_pair(indexBegin(), indexEnd()))
            {
                gapRealigner_.realign(bin_.getBinStart(), bin_.getBinEnd(), index, data_);
            }
            ISAAC_THREAD_CERR << "Realigning gaps done" << std::endl;
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
    const bool keepDuplicates_;
    const alignment::BinMetadata &bin_;
    BamSerializer bamSerializer_;
    std::vector<io::SeFragmentIndex> seIdxFileContent_;
    std::vector<io::RStrandOrShadowFragmentIndex> rIdxFileContent_;
    std::vector<io::FStrandFragmentIndex> fIdxFileContent_;
    PackedFragmentBuffer data_;
    io::FileBufCache<io::FileBufWithReopen> fileBuf_;
    const GapRealignerMode realignGaps_;
    GapRealigner gapRealigner_;
    alignment::BinDataDistribution dataDistribution_;

    void loadData();
    void loadUnalignedData();
    void loadAlignedData();
    const io::FragmentAccessor &loadFragment(std::istream &isData, unsigned long &offset);
    bool isUnalignedBin() const {return bin_.isUnalignedBin();}
    unsigned long getUniqueRecordsCount() const {return isUnalignedBin() ? bin_.getTotalElements() : size();}

    BaseType::iterator indexBegin() {return begin();}
    BaseType::iterator indexEnd() {return end();}

    template<typename InputRecordType>
    void loadIndex(
        const boost::filesystem::path& idxFilePath,
        const unsigned long elements,
        std::vector<InputRecordType> &whereTo)
    {
        ISAAC_THREAD_CERR << "Loading bin index from " << idxFilePath << std::endl;
        const clock_t startLoad = clock();

        std::istream isFIdx(fileBuf_.get(idxFilePath));
        if (!isFIdx) {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open " + idxFilePath.string()));
        }

        whereTo.resize(elements);
        if (!isFIdx.read(reinterpret_cast<char*>(&whereTo.front()), whereTo.size() * sizeof(InputRecordType)))
        {
            BOOST_THROW_EXCEPTION(
                common::IoException(errno,
                                    (boost::format("Failed to read expected %d records from %s") %
                                        whereTo.size() % idxFilePath.string()).str()));
        }
        ISAAC_THREAD_CERR << "Loading bin index done from " << idxFilePath << " in " << (clock() - startLoad) / 1000 << "ms" << std::endl;
    }

};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BIN_SORTER_HH
