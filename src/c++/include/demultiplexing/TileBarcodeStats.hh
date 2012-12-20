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
 ** \file TileBarcodeStats.hh
 **
 ** \brief Statistics collection helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_DEMULTIPLEXING_TILE_BARCODE_STATS_H
#define ISAAC_DEMULTIPLEXING_TILE_BARCODE_STATS_H

namespace isaac
{
namespace demultiplexing
{

struct TileBarcodeStats
{
    unsigned long barcodeCount_;
    unsigned long perfectBarcodeCount_;
    unsigned long oneMismatchBarcodeCount_;

    explicit TileBarcodeStats() :
        barcodeCount_(0), perfectBarcodeCount_(0), oneMismatchBarcodeCount_(0)
    {
    }

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

    const TileBarcodeStats &operator +=(const TileBarcodeStats &right)
    {
        barcodeCount_ += right.barcodeCount_;
        perfectBarcodeCount_ += right.perfectBarcodeCount_;
        oneMismatchBarcodeCount_ += right.oneMismatchBarcodeCount_;
        return *this;
    }

    const TileBarcodeStats operator +(const TileBarcodeStats &right) const
    {
        TileBarcodeStats ret(*this);
        ret += right;
        return ret;
    }

    void finalize()
    {
    }

};

} //namespace demultiplexing
} //namespace isaac

#endif //ISAAC_DEMULTIPLEXING_TILE_BARCODE_STATS_H
