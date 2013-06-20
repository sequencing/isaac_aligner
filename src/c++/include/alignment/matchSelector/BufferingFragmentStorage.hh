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
 ** \file BufferingFragmentStorage.hh
 **
 ** \brief Fragment buffer flushing and output file management.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BUFFERING_FRAGMENT_STORAGE_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_BUFFERING_FRAGMENT_STORAGE_HH

#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/thread.hpp>

#include "alignment/MatchDistribution.hh"
#include "common/Threads.hpp"
#include "io/FileBufCache.hh"

#include "BinIndexMap.hh"
#include "FragmentCollector.hh"
#include "FragmentStorage.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

namespace bfs = boost::filesystem;

class BufferingFragmentStorage: boost::noncopyable, public FragmentStorage
{
public:
    BufferingFragmentStorage(
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
        const unsigned long totalTiles,
        const bool skipEmptyBins);

    virtual void close(alignment::BinMetadataList &binPathList) {binPathList_.swap(binPathList);}

    virtual void add(const BamTemplate &bamTemplate, const unsigned barcodeIdx)
    {
        for (unsigned k = 0; k < bamTemplate.getFragmentCount(); ++k)
        {
            fragmentCollector_.add(bamTemplate, k, barcodeIdx);
        }
    }

    virtual void prepareFlush();
    virtual void flush();
    virtual void resize(const unsigned long clusters)
    {
        fragmentCollector_.resize(clusters);
    }
    virtual void unreserve()
    {
        flushBuffer_.unreserve();
        fragmentCollector_.swapBuffer(flushBuffer_);
        flushBuffer_.unreserve();
    }

private:
    const bool keepUnaligned_;
    const unsigned long maxTileReads_;
    // 0-based number of tile in the order in which they get stored.
    unsigned storedTile_;
    boost::mutex binFlushMutex_;

    const BinIndexMap binIndexMap_;
    common::ThreadVector flushThreads_;

    /// association of a bin index to a path
    alignment::BinMetadataList binPathList_;
    FragmentCollector fragmentCollector_;
    matchSelector::FragmentBuffer flushBuffer_;

    typedef io::FileBufCache<io::FileBufWithReopen > FileBufCache;
    std::vector<FileBufCache> threadDataFileBufCaches_;

    friend std::ostream& operator << (std::ostream& os, const BufferingFragmentStorage &storage);

    void flushBin(
        const unsigned threadNumber,
        const unsigned binNumber, BinMetadata &binMetadata);


    void flushUnmappedBin(
        const unsigned threadNumber,
        const unsigned binNumber, BinMetadata &binMetadata);

    void threadFlushBins(const unsigned threadNumber, unsigned &nextUnflushed);

    /// helper method to initialize the binPathList_
    static alignment::BinMetadataList buildBinPathList(
        const BinIndexMap &binIndexMap,
        const MatchDistribution &matchDistribution,
        const bfs::path &binDirectory,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const unsigned long maxTileClusters,
        const unsigned long totalTiles,
        const bool preSortBins,
        const bool skipEmptyBins);
};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BUFFERING_FRAGMENT_STORAGE_HH
