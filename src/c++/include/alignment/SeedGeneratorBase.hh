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
#include "flowcell/TileMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "io/BclMapper.hh"
#include "oligo/Kmer.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Encapsulates the variables that are shared by all the threads while
 ** loading the seeds.
 **/
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
        const flowcell::ReadMetadataList &readMetadataList,
        const std::vector<SeedMetadata> &seedMetadataList,
        const reference::SortedReferenceXmlList &sortedReferenceXmlList,
        const flowcell::TileMetadataList &tileMetadataList);

    /**
     * \brief returned iterators of the vector point past the last tile for each of the references
     */
    const std::vector<std::vector<Seed>::iterator> &getReferenceSeedBounds() const
    {
        return nextTileSeedBegins_;
    }

protected:
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::ReadMetadataList &readMetadataList_;
    /// count of seeds on each read (seedCounts_[readIndex])
    const std::vector<unsigned> seedCounts_;

    typedef std::vector<std::vector<std::vector<unsigned> > > FragmentCounts;
    FragmentCounts referenceTileReadFragmentCounts_;
    const std::vector<SeedMetadata> seedMetadataOrderedByFirstCycle_;

    /**
     * \brief Geometry: [reference]
     */
    std::vector<std::vector<Seed>::iterator> nextTileSeedBegins_;

    void reset(const flowcell::TileMetadataList &tiles, std::vector<Seed> &seeds,
               const matchFinder::TileClusterInfo &tileClusterBarcode);

    void advanceToNextTile(const flowcell::TileMetadata &currentTile);
    void sortSeeds(
        std::vector<Seed> &seeds,
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
