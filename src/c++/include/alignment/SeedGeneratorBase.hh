/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** BSD 2-Clause License
 **
 ** You should have received a copy of the BSD 2-Clause License
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** \file SeedGeneratorBase.hh
 **
 ** \brief  Helper base class for seed generators
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_SEED_GENERATOR_BASE_HH
#define iSAAC_ALIGNMENT_SEED_GENERATOR_BASE_HH

#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <boost/thread/mutex.hpp>

#include "alignment/Seed.hh"
#include "alignment/SeedMetadata.hh"
#include "common/Threads.hpp"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/TileMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "oligo/Kmer.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Encapsulates the variables that are shared by all the threads while
 ** loading the seeds.
 **/
template <typename KmerT>
class SeedGeneratorBase
{
public:
    /**
     ** \brief constructs an instance with all the required shorthands.
     **
     ** Note: all parameters are kept as references and it is the
     ** responsibility of the caller to ensure appropriate life time for the
     ** referenced variables.
     **/
    SeedGeneratorBase(
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::Layout &flowcellLayout,
        const std::vector<SeedMetadata> &seedMetadataList,
        const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const flowcell::TileMetadataList &tileMetadataList);

    /**
     * \brief returned iterators of the vector point past the last tile for each of the references
     */
    const std::vector<typename std::vector<Seed<KmerT> >::iterator> &getReferenceSeedBounds() const
    {
        return nextTileSeedBegins_;
    }

protected:
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::Layout &flowcellLayout_;
    /// count of seeds on each read (seedCounts_[readIndex])
    const std::vector<unsigned> seedCounts_;

    typedef std::vector<std::vector<std::vector<unsigned> > > FragmentCounts;
    FragmentCounts referenceTileReadFragmentCounts_;
    const std::vector<SeedMetadata> seedMetadataOrderedByFirstCycle_;

    /**
     * \brief Geometry: [reference]
     */
    typename std::vector<typename std::vector<Seed<KmerT> >::iterator> nextTileSeedBegins_;

    void reset(const flowcell::TileMetadataList &tiles, std::vector<Seed<KmerT> > &seeds,
               const matchFinder::TileClusterInfo &tileClusterBarcode);

    void advanceToNextTile(const flowcell::TileMetadata &currentTile);
    void sortSeeds(
        std::vector<Seed<KmerT> > &seeds,
        common::ScoopedMallocBlock  &mallocBlock);

private:
    /// Return the count of seeds for each read
    std::vector<unsigned> getSeedCounts(
        const std::vector<flowcell::ReadMetadata> &readMetadataList,
        const std::vector<SeedMetadata> &seedMetadataList) const;

    void getReferenceTileReadFragmentCounts(
        const flowcell::TileMetadataList &tiles,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const matchFinder::TileClusterInfo &tileClusterBarcode,
        FragmentCounts &fragmentCounts) const;

    static std::vector<SeedMetadata> orderSeedMetadataByFirstCycle(
        const std::vector<SeedMetadata> seedMetadataList);

};


} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SEED_GENERATOR_BASE_HH
