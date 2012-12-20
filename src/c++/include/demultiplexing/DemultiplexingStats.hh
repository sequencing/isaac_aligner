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

#include "demultiplexing/TileBarcodeStats.hh"

namespace isaac
{
namespace demultiplexing
{

class DemultiplexingStats
{
public:
    typedef std::vector<std::pair<Kmer, unsigned long> > UnknownBarcodeHits;
    struct LaneBarcodeStats
    {
        LaneBarcodeStats() {topUnknownBarcodes_.reserve(topUnknownBarcodesMax_);}
        UnknownBarcodeHits topUnknownBarcodes_;
    };

    static const unsigned lanesPerFlowcellMax_ = 8;

private:
    static const unsigned topUnknownBarcodesMax_ = 10;
    static const unsigned totalTilesMax_ = 1000;
    const std::vector<flowcell::BarcodeMetadata> &barcodeMetadataList_;

    std::vector<TileBarcodeStats> tileBarcodeStats_;

    UnknownBarcodeHits topUnknownBarcodes_;

    std::vector<LaneBarcodeStats> laneBarcodeStats_;

public:

    DemultiplexingStats(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList) :
            barcodeMetadataList_(barcodeMetadataList),
            laneBarcodeStats_(flowcellLayoutList.size() * lanesPerFlowcellMax_)
    {
        topUnknownBarcodes_.reserve(topUnknownBarcodesMax_ * 2);
        // TODO: this allows for a situation where every barcode is expected to be found on every tile.
        // Since barcodes are constrained to one lane, plenty of ram can be saved by coming up with an alternative solution.
        const unsigned tileBarcodeStatsCount = barcodeMetadataList.size() * totalTilesMax_;
        ISAAC_THREAD_CERR << "Allocating " << tileBarcodeStatsCount << " demultiplexing tile barcode stats." << std::endl;
        tileBarcodeStats_.resize(tileBarcodeStatsCount);
        ISAAC_THREAD_CERR << "Allocating " << tileBarcodeStatsCount << " demultiplexing tile barcode stats done. Total size is " << tileBarcodeStats_.capacity() * sizeof(TileBarcodeStats) << " bytes."<< std::endl;
    }

    void recordBarcode(const BarcodeId barcodeId)
    {
        tileBarcodeStats_.at(tileBarcodeIndex(barcodeId)).recordBarcode(barcodeId);
    }

    void recordUnknownBarcode(
        const unsigned barcodeIndex,
        const unsigned tile)
    {
        tileBarcodeStats_.at(tileBarcodeIndex(barcodeIndex, tile)).recordUnknownBarcode();
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
        if (topUnknownBarcodesMax_ == topUnknownBarcodes_.size())
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
     * \param lane 1-based lane number
     */
    void finalizeUnknownBarcodeHits(const unsigned flowcellIndex, const unsigned lane)
    {
        ISAAC_ASSERT_MSG(lane <= lanesPerFlowcellMax_, "Unexpected lane number");
        LaneBarcodeStats &laneStats = laneBarcodeStats_.at(flowcellIndex * lanesPerFlowcellMax_ + lane - 1);
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
        const unsigned topUnknownBarcodesMax = topUnknownBarcodesMax_;
        topUnknownBarcodes_.resize(std::min<unsigned>(topUnknownBarcodesMax,
                                                      std::distance(topUnknownBarcodes_.begin(), first)));
        laneStats.topUnknownBarcodes_.clear();
        laneStats.topUnknownBarcodes_.insert(laneStats.topUnknownBarcodes_.end(), topUnknownBarcodes_.begin(), topUnknownBarcodes_.end());
        topUnknownBarcodes_.clear();
    }
/*
    DemultiplexingStats &operator +=(const DemultiplexingStats &right)
    {
        ISAAC_ASSERT_MSG(right.barcodeMetadataList_.size() == barcodeMetadataList_.size(), "dimensions must match");
        ISAAC_ASSERT_MSG(right.tileBarcodeStats_.size() == tileBarcodeStats_.size(), "size must match");
        std::transform(tileBarcodeStats_.begin(), tileBarcodeStats_.end(),
                       right.tileBarcodeStats_.begin(), tileBarcodeStats_.begin(), std::plus<TileBarcodeStats>());
        return *this;
    }

    const DemultiplexingStats operator +(const DemultiplexingStats &right) const
    {
        DemultiplexingStats ret(*this);
        ret += right;
        return ret;
    }

    DemultiplexingStats & operator =(const DemultiplexingStats &that) {
        ISAAC_ASSERT_MSG(that.barcodeMetadataList_.size() == barcodeMetadataList_.size(), "dimensions must match");
        ISAAC_ASSERT_MSG(that.tileBarcodeStats_.size() == tileBarcodeStats_.size(), "size must match");
        tileBarcodeStats_ = that.tileBarcodeStats_;
        return *this;
    }*/

    const TileBarcodeStats &getTileBarcodeStat(
        const flowcell::BarcodeMetadata& barcode,
        const flowcell::TileMetadata& tile) const
    {
        return tileBarcodeStats_.at(tileBarcodeIndex(barcode.getIndex(), tile.getIndex()));
    }

    /**
     * \param lane 1-based lane number
     */
    const LaneBarcodeStats &getFlowcellLaneStat(
        const flowcell::Layout &flowcell,
        const unsigned lane) const
    {
        return laneBarcodeStats_.at(flowcell.getIndex() * lanesPerFlowcellMax_ + lane - 1);
    }

private:

    unsigned tileBarcodeIndex(
        const demultiplexing::BarcodeId& barcode) const
    {
        return tileBarcodeIndex(barcode.getBarcode(), barcode.getTile());
    }

    unsigned tileBarcodeIndex(
        const flowcell::BarcodeMetadata& barcode,
        const flowcell::TileMetadata& tile) const
    {
        return tileBarcodeIndex(barcode.getIndex(), tile.getIndex());
    }

    unsigned tileBarcodeIndex(const unsigned barcode, const unsigned tile) const
    {
        ISAAC_ASSERT_MSG(totalTilesMax_ > tile, "Maximum number of tiles exceeded in DemultiplexingStats");
        return tile * barcodeMetadataList_.size() + barcode;
    }
};

} //namespace demultiplexing
} //namespace isaac

#endif //ISAAC_DEMULTIPLEXING_DEMULTIPLEXING_STATS_H
