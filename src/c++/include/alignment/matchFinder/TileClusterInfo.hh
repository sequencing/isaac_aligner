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
 ** \file TileClusterBarcode.hh
 **
 ** \brief Holds the information about cluster barcode mapping.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_TILE_CLUSTER_INFO_HH
#define iSAAC_DEMULTIPLEXING_TILE_CLUSTER_INFO_HH

#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <boost/thread/mutex.hpp>

#include "alignment/Seed.hh"
#include "alignment/SeedMetadata.hh"
#include "common/Threads.hpp"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "oligo/Kmer.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace alignment
{
namespace matchFinder
{

/**
 * \brief Contains the cluster information required across multiple passes of match finding
 *
 * At the moment the following data is required:
 *  1. cluster barcode index
 *  2. per-read indicator of whether the further match finding is needed
 *  3. last bit of each byte is currently unused
 *
 * The read information is split across between two different memory bytes to allow for asynchronous
 * updates from different threads.
 * Layout:
 *
 *  +-+-------------+-+-------------+
 *  |0|1|2 3 4 5 6 7|0|1 2 3 4 5 6 7|
 *  +-+-------------+-+-------------+
 *  |    byte1_     |    byte2_     |
 *  +-+-----------+-+-+-----------+-+
 *  |r|barcode    |?|r|barcode    |?|
 *  |1|idx part 1 |?|2|idx part 2 |?|
 *  +-+-----------+-+-+-----------+-+
 */
class ClusterInfo
{
public:
    // width in bits for each field
    static const unsigned R1FOUND_WIDTH   = 1;
    static const unsigned BARCODE1_WIDTH  = 6;
    static const unsigned R2FOUND_WIDTH   = 1;
    static const unsigned BARCODE2_WIDTH  = 6;
    // masks for the values in each field
    static const unsigned char R1FOUND_MASK = ~(~0U<<R1FOUND_WIDTH);
    static const unsigned char BARCODE1_MASK = ~(~0U<<(BARCODE1_WIDTH)) << R1FOUND_WIDTH;
    static const unsigned char R2FOUND_MASK = ~(~0U<<R2FOUND_WIDTH);
    static const unsigned char BARCODE2_MASK = ~(~0U<<(BARCODE2_WIDTH)) << R2FOUND_WIDTH;
    // shifts in bits for each field
    static const unsigned R1FOUND_SHIFT = 0;
    static const unsigned BARCODE1_SHIFT = R1FOUND_SHIFT + R1FOUND_WIDTH;
    static const unsigned R2FOUND_SHIFT = 0;
    static const unsigned BARCODE2_SHIFT = R2FOUND_SHIFT + R2FOUND_WIDTH;

    static const unsigned MAX_BARCODE_VALUE = (BARCODE1_MASK >> BARCODE1_SHIFT) |
                                             ((BARCODE2_MASK >> BARCODE2_SHIFT) << BARCODE1_WIDTH);
    ClusterInfo() : byte1_(-1), byte2_(-1)
    {
        // ensure initially barcode is set to someting we can treat as uninitialized.
        unmarkComplete();
    }
    ClusterInfo(const bool markComplete) :
        byte1_(markComplete ? R1FOUND_MASK : 0), byte2_(markComplete ? R2FOUND_MASK : 0){}

    unsigned getBarcodeIndex() const
    {
        return ((byte1_ & BARCODE1_MASK) >> BARCODE1_SHIFT) |
              (((byte2_ & BARCODE2_MASK) >> BARCODE2_SHIFT) << BARCODE1_WIDTH);
    }

    bool isBarcodeSet() const
    {
        return MAX_BARCODE_VALUE != getBarcodeIndex();
    }

    void setBarcodeIndex(const unsigned barcodeIndex)
    {
        ISAAC_ASSERT_MSG(barcodeIndex < MAX_BARCODE_VALUE, "Barcode does not fit in the allowed bit range");
        byte1_ = (byte1_ & R1FOUND_MASK) | ((barcodeIndex << BARCODE1_SHIFT) & (BARCODE1_MASK));
        byte2_ = (byte2_ & R2FOUND_MASK) | (((barcodeIndex >> BARCODE1_WIDTH) << BARCODE1_SHIFT) & BARCODE2_MASK);
    }

    bool isReadComplete(const unsigned readIndex) const
    {
        if (0 == readIndex)
        {
            return byte1_ & R1FOUND_MASK;
        }
        else
        {
            return byte2_ & R2FOUND_MASK;
        }
    }

    void markReadComplete(const unsigned readIndex)
    {
        if (0 == readIndex)
        {
            byte1_ |= R1FOUND_MASK;
        }
        else
        {
            byte2_ |= R2FOUND_MASK;
        }
    }

    void unmarkComplete()
    {
        byte1_ &= ~R1FOUND_MASK;
        byte2_ &= ~R2FOUND_MASK;
    }

private:
    unsigned char byte1_;
    unsigned char byte2_;
};

inline std::ostream& operator << (std::ostream &os, const ClusterInfo &cluster)
{
    return os << "ClusterInfo(" << cluster.getBarcodeIndex() <<
        ", " << cluster.isReadComplete(0) <<
        ", " << cluster.isReadComplete(1) <<
        ")";
}

/**
 * \brief Geometry: [tileIndex][clusterIndex].
 *
 * For each tile contains mapping between the tile cluster id and the index
 * of the barcode found for the cluster.
 */
struct TileClusterInfo : std::vector<std::vector<ClusterInfo> >
{
    TileClusterInfo(const flowcell::TileMetadataList &unprocessedTileMetadataList,
                    const std::vector<size_t> &clusterIdFilter):
        std::vector<std::vector<ClusterInfo> >(unprocessedTileMetadataList.back().getIndex() + 1)
    {
        if (clusterIdFilter.size())
        {
            // if filter is given, create all as complete and then umark the listed clusters
            BOOST_FOREACH(const flowcell::TileMetadata &unprocessedTile, unprocessedTileMetadataList)
            {
                std::vector<ClusterInfo> &tile = at(unprocessedTile.getIndex());
                tile.resize(unprocessedTile.getClusterCount(), ClusterInfo(true));
                BOOST_FOREACH(const size_t cluster, clusterIdFilter)
                {
                    tile.at(cluster).unmarkComplete();
                }
            }
        }
        else
        {
            // allocated space for only the tiles we'll be working on
            BOOST_FOREACH(const flowcell::TileMetadata &unprocessedTile, unprocessedTileMetadataList)
            {
                std::vector<ClusterInfo> &tile = at(unprocessedTile.getIndex());
                tile.resize(unprocessedTile.getClusterCount());
            }
        }
    }
    unsigned getBarcodeIndex(const unsigned tileIndex, const unsigned clusterIndex) const
    {
        return at(tileIndex).at(clusterIndex).getBarcodeIndex();
    }

    void setBarcodeIndex(const unsigned tileIndex, const unsigned clusterIndex, const unsigned barcodeIndex)
    {
        at(tileIndex).at(clusterIndex).setBarcodeIndex(barcodeIndex);
    }

    bool isReadComplete(const unsigned tileIndex, const unsigned clusterIndex, const unsigned readIndex) const
    {
        return at(tileIndex).at(clusterIndex).isReadComplete(readIndex);
    }

    void markReadComplete(const unsigned tileIndex, const unsigned clusterIndex, const unsigned readIndex)
    {
        at(tileIndex).at(clusterIndex).markReadComplete(readIndex);
    }
};

} // namespace matchFinder
} // namespace demultiplexing
} // namespace isaac

#endif // #ifndef iSAAC_DEMULTIPLEXING_TILE_CLUSTER_INFO_HH
