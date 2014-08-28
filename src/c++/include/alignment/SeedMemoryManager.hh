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
 ** \file SeedMemoryManager.hh
 **
 ** \brief Deals with predicting and allocating memory for seeds.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_SEED_MEMORY_MANAGER_HH
#define iSAAC_ALIGNMENT_SEED_MEMORY_MANAGER_HH

#include <vector>
#include <boost/noncopyable.hpp>

#include "alignment/Seed.hh"
#include "alignment/SeedMetadata.hh"
#include "common/Threads.hpp"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "flowcell/ReadMetadata.hh"

namespace isaac
{
namespace alignment
{

template <typename KmerT>
class SeedMemoryManager: boost::noncopyable
{
    typedef alignment::Seed<KmerT> SeedT;

public:
    typedef flowcell::TileMetadataList TileMetadataList;
    typedef flowcell::ReadMetadataList ReadMetadataList;

    SeedMemoryManager(
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const ReadMetadataList &readMetadataList,
        const SeedMetadataList &seedMetadataList,
        const flowcell::TileMetadataList &unprocessedTileMetadataList);

    bool selectTiles(TileMetadataList &unprocessedPool,
                     const matchFinder::TileClusterInfo &fragmentsToSkip,
                     const unsigned maxTilesAtATime,
                     const unsigned maxSavers,
                     TileMetadataList &selectedTiles);

    void allocate(const TileMetadataList &tiles, std::vector<SeedT> &seeds) const;

private:
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const ReadMetadataList &readMetadataList_;
    const SeedMetadataList &seedMetadataList_;

    // notFoundMatchesCount_[readIndex][tileIndex] : count of fragmentsToSkip[readIndex][tileIndex][clusterId] != true
    std::vector<std::vector<unsigned > > notFoundMatchesCount_;

    unsigned long getTotalSeedCount(const TileMetadataList &tiles) const;
    const std::vector<std::vector<unsigned > > getNotFoundMatchesCount(
        const flowcell::TileMetadataList &unprocessedTiles,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const ReadMetadataList &readMetadataList,
        const matchFinder::TileClusterInfo &foundMatches) const;

    bool seeIfFits(const TileMetadataList &tiles) const;
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SEED_MEMORY_MANAGER_HH
