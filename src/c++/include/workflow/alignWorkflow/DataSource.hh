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
 ** \file DataSource.hh
 **
 ** \brief Abstraction of data source
 **
 ** \author Roman Petrovski
 **/


#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_DATA_SOURCE_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_DATA_SOURCE_HH

#include "alignment/BclClusters.hh"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

template <typename KmerT>
struct SeedSource : boost::noncopyable
{
    typedef alignment::Seed<KmerT> SeedT;
    typedef typename std::vector<SeedT>::iterator SeedIterator;
    /**
     * \brief Returns set of tiles that can be processed together.
     *
     * \return If returned set is empty, there is nothing left to process
     */
    virtual flowcell::TileMetadataList discoverTiles() = 0;

    /**
     * \brief Initializes internal buffers, remembers information needed for seed generation
     */
    virtual void initBuffers(
        flowcell::TileMetadataList &unprocessedTiles,
        const alignment::SeedMetadataList &seedMetadataList,
        common::ThreadVector &threads) = 0;

    /**
     * \brief Generates seeds for tiles based on seedMetadataList supplied to initBuffers
     */
    virtual void generateSeeds(
        const flowcell::TileMetadataList &tiles,
        const alignment::matchFinder::TileClusterInfo &tileClusterBarcode,
        std::vector<SeedT> &seeds,
        common::ScoopedMallocBlock  &mallocBlock) = 0;

    /**
     * \brief Returns the list of iterators denoting the boundaries of seeds that must be aligned
     *        against the corresponding reference genome. The exchange of reference metadata between
     *        client code and implementation is implementation-specific.
     */
    virtual const std::vector<SeedIterator> &getReferenceSeedBounds() const = 0;

    virtual ~SeedSource(){}
};

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac


#endif //#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_DATA_SOURCE_HH
