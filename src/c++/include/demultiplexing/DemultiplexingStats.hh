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
 ** \file DemultiplexingStats.hh
 **
 ** \brief BarcodeResolver statistics helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_DEMULTIPLEXING_DEMULTIPLEXING_STATS_H
#define ISAAC_DEMULTIPLEXING_DEMULTIPLEXING_STATS_H

#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace demultiplexing
{

static const unsigned TOP_UNKNOWN_BARCODES_MAX = 10;
typedef std::vector<std::pair<Kmer, unsigned long> > UnknownBarcodeHits;
struct LaneBarcodeStats
{
    UnknownBarcodeHits topUnknownBarcodes_;
    unsigned long barcodeCount_;
    unsigned long perfectBarcodeCount_;
    unsigned long oneMismatchBarcodeCount_;

    LaneBarcodeStats():barcodeCount_(0), perfectBarcodeCount_(0), oneMismatchBarcodeCount_(0)
        {topUnknownBarcodes_.reserve(TOP_UNKNOWN_BARCODES_MAX);}

    void recordBarcode(const BarcodeId &barcodeId)
    {
        ++barcodeCount_;
        perfectBarcodeCount_ += (0 == barcodeId.getMismatches());
        oneMismatchBarcodeCount_ += (1 == barcodeId.getMismatches());
    }

    void recordUnknownBarcode()
    {
        ++barcodeCount_;
    }

    const LaneBarcodeStats &operator +=(const LaneBarcodeStats &right)
    {
        barcodeCount_ += right.barcodeCount_;
        perfectBarcodeCount_ += right.perfectBarcodeCount_;
        oneMismatchBarcodeCount_ += right.oneMismatchBarcodeCount_;
        return *this;
    }

};

class DemultiplexingStats
{
private:
    static const unsigned TOTAL_TILES_MAX = 1000;
    const std::vector<flowcell::BarcodeMetadata> &barcodeMetadataList_;

    UnknownBarcodeHits topUnknownBarcodes_;

    std::vector<LaneBarcodeStats> laneBarcodeStats_;

public:

    DemultiplexingStats(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList) :
            barcodeMetadataList_(barcodeMetadataList),
            laneBarcodeStats_(barcodeMetadataList_.size())
    {
        topUnknownBarcodes_.reserve(TOP_UNKNOWN_BARCODES_MAX * 2);
    }

    void recordBarcode(const BarcodeId barcodeId)
    {
        laneBarcodeStats_.at(laneBarcodeIndex(barcodeId)).recordBarcode(barcodeId);
    }

    void recordUnknownBarcode(
        const unsigned barcodeIndex,
        const unsigned tile)
    {
        laneBarcodeStats_.at(laneBarcodeIndex(barcodeIndex)).recordUnknownBarcode();
    }

    static bool orderBySequence(
        const std::pair<Kmer, unsigned long> &left,
        const std::pair<Kmer, unsigned long> &right)
    {
        return left.first < right.first;
    }

    static bool orderByHits(
        const std::pair<Kmer, unsigned long> &left,
        const std::pair<Kmer, unsigned long> &right)
    {
        // we want the Greatest Hits on top!
        return left.second > right.second;
    }

    void recordUnknownBarcodeHits(const Kmer sequence, unsigned long hits)
    {
        UnknownBarcodeHits::iterator insertIt =
            std::lower_bound(topUnknownBarcodes_.begin(), topUnknownBarcodes_.end(),
                             std::make_pair(0UL, hits), orderByHits);
        if (TOP_UNKNOWN_BARCODES_MAX == topUnknownBarcodes_.size())
        {
            if (topUnknownBarcodes_.end() != insertIt)
            {
                topUnknownBarcodes_.pop_back();
            }
            else
            {
                return;
            }
        }
        topUnknownBarcodes_.insert(insertIt, std::make_pair(sequence, hits));
    }

    /**
     * \brief merge with already accumulated for lane and extract the  top popular ones
     *
     * \param barcode index of the 'unknown' barcode of that lane
     */
    void finalizeUnknownBarcodeHits(const unsigned barcodeIndex)
    {
        ISAAC_ASSERT_MSG(barcodeMetadataList_.at(barcodeIndex).isUnknown(), "Barcode index does not designate lane unknown barcode" << barcodeMetadataList_.at(barcodeIndex));
        LaneBarcodeStats &laneStats = laneBarcodeStats_.at(laneBarcodeIndex(barcodeIndex));
        topUnknownBarcodes_.insert(topUnknownBarcodes_.end(), laneStats.topUnknownBarcodes_.begin(), laneStats.topUnknownBarcodes_.end());
        std::sort(topUnknownBarcodes_.begin(), topUnknownBarcodes_.end(), orderBySequence);
        UnknownBarcodeHits::iterator first = topUnknownBarcodes_.begin();
        UnknownBarcodeHits::iterator second = topUnknownBarcodes_.begin();
        while (second != topUnknownBarcodes_.end())
        {
            if (++second != topUnknownBarcodes_.end())
            {
                if (first->first == second->first)
                {
                    first->second += second->second;
                }
                else
                {
                    *(++first) = *second;
                }
            }
            else
            {
                // make sure first points past the last valid element
                ++first;
                break;
            }
        }
        std::sort(topUnknownBarcodes_.begin(), first, orderByHits);
        // else gcc fails at link time with -O0
        const unsigned topUnknownBarcodesMax = TOP_UNKNOWN_BARCODES_MAX;
        topUnknownBarcodes_.resize(std::min<unsigned>(topUnknownBarcodesMax,
                                                      std::distance(topUnknownBarcodes_.begin(), first)));
        laneStats.topUnknownBarcodes_.clear();
        laneStats.topUnknownBarcodes_.insert(laneStats.topUnknownBarcodes_.end(), topUnknownBarcodes_.begin(), topUnknownBarcodes_.end());
        topUnknownBarcodes_.clear();
    }

    const LaneBarcodeStats &getLaneBarcodeStat(
        const flowcell::BarcodeMetadata& barcode) const
    {
        return laneBarcodeStats_.at(laneBarcodeIndex(barcode.getIndex()));
    }

    const LaneBarcodeStats &getLaneUnknwonBarcodeStat(
        const unsigned barcodeIndex) const
    {
        ISAAC_ASSERT_MSG(barcodeMetadataList_.at(barcodeIndex).isUnknown(), "Barcode index does not designate lane unknown barcode" << barcodeMetadataList_.at(barcodeIndex));
        return laneBarcodeStats_.at(laneBarcodeIndex(barcodeIndex));
    }

private:

    unsigned laneBarcodeIndex(
        const demultiplexing::BarcodeId& barcode) const
    {
        return laneBarcodeIndex(barcode.getBarcode());
    }

    unsigned laneBarcodeIndex(const unsigned barcode) const
    {
        return barcode;
    }
};

} //namespace demultiplexing
} //namespace isaac

#endif //ISAAC_DEMULTIPLEXING_DEMULTIPLEXING_STATS_H
