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
 ** \file SeedLoader.hh
 **
 ** \brief Loads the seeds from the BCL files in parallel.
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_SEED_LOADER_HH
#define iSAAC_ALIGNMENT_SEED_LOADER_HH

#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <boost/thread/mutex.hpp>

#include "alignment/Seed.hh"
#include "alignment/SeedGeneratorBase.hh"
#include "alignment/SeedMetadata.hh"
#include "common/Threads.hpp"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "oligo/Kmer.hh"
#include "rta/BclMapper.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Encapsulates the variables that are shared by all the threads while
 ** loading the seeds.
 **/
template <typename ReaderT, typename KmerT>
class ParallelSeedLoader : private SeedGeneratorBase<KmerT>
{
    typedef SeedGeneratorBase<KmerT> BaseT;
public:
    /**
     ** \brief constructs an instance with all the required shorthands.
     **
     ** Note: all parameters are kept as references and it is the
     ** responsibility of the caller to ensure appropriate life time for the
     ** referenced variables.
     **/
    ParallelSeedLoader(
        const bool ignoreMissingBcls,
        common::ThreadVector &threads,
        boost::ptr_vector<rta::SingleCycleBclMapper<ReaderT> > &threadBclMappers,
        const unsigned inputLoadersMax,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::Layout &flowcellLayout,
        const std::vector<SeedMetadata> &seedMetadataList,
        const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const flowcell::TileMetadataList &tileMetadataList);

    void loadSeeds(const flowcell::TileMetadataList &tiles,
                   const matchFinder::TileClusterInfo &tileClusterBarcode,
                   std::vector<Seed<KmerT> > &seeds,
                   common::ScoopedMallocBlock  &mallocBlock);

    /**
     * \brief returned iterators of the vector point past the last tile for each of the references
     */
    const std::vector<typename std::vector<Seed<KmerT> >::iterator> &getReferenceSeedBounds() const
    {
        return BaseT::nextTileSeedBegins_;
    }

private:
    // The mutex used to acquire the next tile and the destination of the seeds
    boost::mutex mutex_;
    const unsigned inputLoadersMax_;

    const std::vector<unsigned> seedCycles_;

    /**
     * \brief Geometry: [thread][reference]
     */
    std::vector<std::vector<typename std::vector<Seed<KmerT> >::iterator> > threadDestinations_;
    std::vector<std::vector<typename std::vector<Seed<KmerT> >::iterator> > threadCycleDestinations_;

    common::ThreadVector &threads_;
    boost::ptr_vector<rta::SingleCycleBclMapper<ReaderT> > &threadBclMappers_;

    void load(
        const matchFinder::TileClusterInfo &tileClusterBarcode,
        std::vector<flowcell::TileMetadata>::const_iterator &nextTile,
        const std::vector<flowcell::TileMetadata>::const_iterator tilesEnd,
        const unsigned threadNumber);

    void loadTileCycle(
        const matchFinder::TileClusterInfo &tileClusterBarcode,
        rta::SingleCycleBclMapper<ReaderT> &bclMapper,
        std::vector<typename std::vector<Seed<KmerT> >::iterator> &destinationBegins,
        const flowcell::TileMetadata &tile,
        const unsigned cycle,
        std::vector<SeedMetadata>::const_iterator cycleSeedsBegin,
        const std::vector<SeedMetadata>::const_iterator cycleSeedsEnd);

    std::vector<unsigned> discoverSeedCycles() const;
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SEED_LOADER_HH
