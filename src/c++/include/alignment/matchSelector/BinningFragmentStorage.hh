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
 ** \file BinningFragmentStorage.hh
 **
 ** \brief Stores fragments in bin files without buffering.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BINNING_FRAGMENT_STORAGE_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_BINNING_FRAGMENT_STORAGE_HH

#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/thread.hpp>

#include "alignment/MatchDistribution.hh"
#include "common/Threads.hpp"
#include "io/FileBufCache.hh"

#include "BinIndexMap.hh"
#include "FragmentStorage.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

namespace bfs = boost::filesystem;

class BinningFragmentStorage: boost::noncopyable, public FragmentStorage
{
public:
    BinningFragmentStorage(
        const bool keepUnaligned,
        const bool preSortBins,
        const unsigned maxSavers,
        const unsigned threadBuffers,
        const MatchDistribution &matchDistribution,
        const unsigned long outputBinSize,
        const bfs::path &binDirectory,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const unsigned long maxTileClusters,
        const unsigned long totalTiles);

    virtual void close(alignment::BinMetadataList &binPathList) {binPathList_.swap(binPathList);}

    virtual void add(const BamTemplate &bamTemplate, const unsigned barcodeIdx);

    virtual void prepareFlush()
    {
    }
    virtual void flush()
    {
    }
    virtual void resize(const unsigned long clusters)
    {
    }
    virtual void unreserve()
    {
    }

private:
    static const unsigned READS_MAX = 2;
    const bool keepUnaligned_;
    const unsigned long maxTileReads_;

    const BinIndexMap binIndexMap_;

    /// association of a bin index to a path
    alignment::BinMetadataList binPathList_;
    boost::array<boost::mutex, 8> binMutex_;
    boost::ptr_vector<std::ofstream > binFiles_;

    friend std::ostream& operator << (std::ostream& os, const BinningFragmentStorage &storage);

    /// helper method to initialize the binPathList_
    static alignment::BinMetadataList buildBinPathList(
        const BinIndexMap &binIndexMap,
        const unsigned long outputBinSize,
        const bfs::path &binDirectory,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const unsigned long maxTileClusters,
        const unsigned long totalTiles,
        const bool preSortBins);

    template <typename BufferT>
    void storeFragment(
        const BufferT &buffer,
        const unsigned storageBin);

    template <typename BufferT>
    unsigned packFragment(
        const alignment::BamTemplate &bamTemplate,
        const unsigned fragmentIndex,
        const unsigned barcodeIdx,
        BufferT &buffer);

    void storePaired(const BamTemplate &bamTemplate, const unsigned barcodeIdx);
    void storeSingle(const BamTemplate &bamTemplate, const unsigned barcodeIdx);
};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BINNING_FRAGMENT_STORAGE_HH
