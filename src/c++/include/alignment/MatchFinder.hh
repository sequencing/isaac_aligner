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
 ** \file MatchFinder.hh
 **
 ** \brief Find all matches between seeds and the complete reference.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_FINDER_HH
#define iSAAC_ALIGNMENT_MATCH_FINDER_HH

#include <string>
#include <vector>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/noncopyable.hpp>
#include <boost/thread/mutex.hpp>

#include "alignment/SeedMetadata.hh"
#include "alignment/SeedId.hh"
#include "alignment/Seed.hh"
#include "alignment/MatchTally.hh"
#include "alignment/MatchDistribution.hh"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "common/Threads.hpp"
#include "flowcell/TileMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "io/FileBufCache.hh"
#include "io/MatchWriter.hh"
#include "oligo/Kmer.hh"
#include "reference/ReferenceKmer.hh"
#include "statistics/MatchFinderTileStats.hh"

namespace isaac
{
namespace alignment
{

namespace bfs = boost::filesystem;

template <typename KmerT>
class MatchFinder: boost::noncopyable
{
public:
    typedef Seed<KmerT> SeedT;
    typedef reference::ReferenceKmer<KmerT> ReferenceKmer;

    MatchFinder(
        const reference::SortedReferenceMetadataList &sortedReferenceList,
        const bfs::path & tempDirectory,
        const flowcell::TileMetadataList &tiles,
        const flowcell::ReadMetadataList &readMetadataList,
        const SeedMetadataList &seedMetadataList,
        const unsigned iteration,
        const bool ignoreNeighbors,
        const bool ignoreRepeats,
        const unsigned repeatThreshold,
        const unsigned neighborhoodSizeThreshold,
        MatchTally &matchTally,
        matchFinder::TileClusterInfo &foundMatches,
        common::ThreadVector &threads,
        const unsigned coresMax,
        const unsigned tempSaversMax,
        const unsigned unavailableFileHandles);

    void setTiles(const flowcell::TileMetadataList &tiles);

    /**
     ** \brief Find all the matches for the given list of seeds
     **
     ** The matches can have up to two mismatches. Each N counts as a mismatch.
     ** All the matches are stored in temporary files. The temporary files are
     ** split into several bins. It is the responsibility of the MatchWriter to
     ** decide of the binning.
     **
     ** \param[in] seedIndexList The list of indexes for the seeds to actually
     ** use (as indexed in the SeedMetadataList initially provided to the
     ** MatchFinder constructor)
     **
     ** \param[in] iteration The unique identifier for the current iteration of
     ** match finding (typically the iteration number) used to identify the
     ** match files.
     **
     ** \param repeatThreshold the number of repeats above which the matches are
     ** discarded
     **/
    const std::vector<MatchDistribution> & findMatches(
        const typename std::vector<SeedT>::iterator seedsBegin,
        const std::vector<typename std::vector<SeedT>::iterator> &referenceSeedBounds,
        const bool findNeighbors,
        const bool finalPass);

    unsigned getMaxTileCount() {return maxTilesAtATime_;}
private:
    const std::vector<size_t> clusterIdList_;

    struct KmerSourceMetadata
    {
        KmerSourceMetadata(unsigned referenceIndex,
                           unsigned maskWidth,
                           unsigned mask,
                           boost::filesystem::path maskFilePath):
                               referenceIndex_(referenceIndex), maskWidth_(maskWidth),
                               mask_(mask), maskFilePath_(maskFilePath){}
        unsigned referenceIndex_;
        unsigned maskWidth_;
        unsigned mask_;
        boost::filesystem::path maskFilePath_;

        size_t getPathSize() const {return maskFilePath_.string().size();}
    };
    typedef std::vector<KmerSourceMetadata> KmerSourceMetadataList;
    const KmerSourceMetadataList kmerSourceMetadataList_;

    /**
     * [reference][reference contig index]
     */
    const std::vector<std::vector<unsigned> > referenceContigKaryotypes_;

    const SeedMetadataList &seedMetadataList_;

//    const std::map<std::string, std::vector<boost::filesystem::path> > maskFiles_;
    /// the current iteration
    const unsigned iteration_;
    const bool ignoreNeighbors_;
    // if set, repeat matches will not be reported to match selector, therefore it will try
    // to place fragmets into dodgy locations while boosting the % aligned
    const bool ignoreRepeats_;
    /// the current repeat-threshold to use (number of repeats above which the matches are discarded)
    const unsigned repeatThreshold_;
    const unsigned neighborhoodSizeThreshold_;

    /// foundMatches_[readIndex][tileIndex][clusterId]
    matchFinder::TileClusterInfo &foundExactMatchesOnly_;
    common::ThreadVector &threads_;
    const unsigned threadsMax_;
    const unsigned maxTilesAtATime_;

    std::vector<std::vector<ReferenceKmer> > threadRepeatLists_;
    std::vector<std::vector<ReferenceKmer> > threadNeighborsLists_;
    typedef std::vector<MatchDistribution> ThreadMatchDistributions;
    ThreadMatchDistributions threadMatchDistributions_;

    io::TileMatchWriter matchWriter_;
    std::vector<io::FileBufCache<io::FileBufWithReopen> > threadReferenceFileBuffers_;


    mutable boost::mutex mutex_;

    /// top level component to find all the matches for the currently loaded seeds
    const std::vector<MatchDistribution> & match(
        typename std::vector<SeedT>::const_iterator seedsBegin,
        const std::vector<typename std::vector<SeedT>::iterator> &referenceSeedBounds,
        const bool findNeighbors,
        const bool finalPass);
    /// management of parallelism between concurrent matchMask operations
    void matchMaskParallel(
        typename std::vector<SeedT>::const_iterator &seedsBegin,
        const std::vector<typename std::vector<SeedT>::iterator> &referenceSeedBounds,
        const bool findNeighbors,
        const bool finalPass,
        typename KmerSourceMetadataList::const_iterator &kmerSourceIterator,
        const unsigned threadNumber);

    std::pair<typename std::vector<SeedT>::const_iterator, typename std::vector<SeedT>::const_iterator>
        skipToTheNextMask(
        const typename std::vector<SeedT>::const_iterator currentBegin,
        const typename std::vector<SeedT>::const_iterator seedsEnd,
        const KmerT currentMask,
        const unsigned maskWidth,
        const bool storeNSeedNoMatches);

    const std::vector<KmerSourceMetadata> getMaskFilesList(
        const reference::SortedReferenceMetadataList &sortedReferenceList) const;

    unsigned verifyMaxTileCount(
        const unsigned unavailableFileHandlesCount,
        const unsigned maxSavers,
        const unsigned tilesCount) const;

};

} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_FINDER_HH
