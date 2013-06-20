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
 ** \file Barcode.hh
 **
 ** Barcode identification and manipulation
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_BARCODE_HH
#define iSAAC_DEMULTIPLEXING_BARCODE_HH

#include <boost/format.hpp>
#include <boost/mpl/equal_to.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "oligo/Nucleotides.hh"

namespace isaac
{
namespace demultiplexing
{

/**
 * A - 0
 * C - 1
 * G - 2
 * T - 3
 * N - 4
 *
 * \sa isaac::oligo::getBase
 */
typedef unsigned long Kmer;
static const unsigned int BITS_PER_BASE = 3;
static const Kmer kmerMask_ = ~(Kmer(0)) >> (sizeof(Kmer) * 8 - BITS_PER_BASE);
static const unsigned int MAX_BARCODE_LENGTH = (sizeof(Kmer) * 8) / BITS_PER_BASE;

inline std::string bases(Kmer kmer, unsigned kmerLength)
{
    ISAAC_ASSERT_MSG(MAX_BARCODE_LENGTH >= kmerLength, "kmerLength must be no greater than total number o bases a barcode can fit");
    return oligo::bases<BITS_PER_BASE>(kmer, kmerLength);
}


inline std::ostream & operator <<(std::ostream &os, const oligo::Bases<BITS_PER_BASE, Kmer> &bases)
{
    return oligo::printBases(os, bases);
}


/**
 ** \brief Structured unique identifier of a barcode.
 **
 ** From LSB to MSB:
 **   - cluster:    31 (2,147,483,648)
 **   - barcode:    12 (        4,096)
 **   - tile   :    12 (        4,096)
 **   - mismatches   2 (            4)
 **
 **
 ** The order of the fields is important as it will define the natural order of
 ** the barcode ids when sorted.
 **
 ** \note The implementation choice at the moment is to use a tile number that
 ** is the sequential index (0-based) of the tile in the input data set. This
 ** solution has been selected instead of keeping the actual tile and lane, to
 ** make easier (more flexible) to retrieve the metadata and components
 ** associated to each tile throughout the application.
 **
 ** \sa isaac::alignment::SeedId
 **/
class BarcodeId
{
public:
    // width in bits for each field
    static const unsigned MISMATCHES_WIDTH = 2;
    static const unsigned CLUSTER_WIDTH = 31;
    static const unsigned BARCODE_WIDTH = 12;
    static const unsigned TILE_WIDTH = 12;
    // masks for the values in each field
    static const unsigned long MISMATCHES_MASK = ~(~0UL<<MISMATCHES_WIDTH);
    static const unsigned long CLUSTER_MASK = ~(~0UL<<CLUSTER_WIDTH);
    static const unsigned long BARCODE_MASK = ~(~0UL<<BARCODE_WIDTH);
    static const unsigned long TILE_MASK = ~(~0UL<<TILE_WIDTH);
    // shifts in bits for each field
    static const unsigned MISMATCHES_SHIFT = 0;
    static const unsigned CLUSTER_SHIFT = MISMATCHES_SHIFT + MISMATCHES_WIDTH;
    static const unsigned BARCODE_SHIFT = CLUSTER_SHIFT + CLUSTER_WIDTH;
    static const unsigned TILE_SHIFT = BARCODE_SHIFT + BARCODE_WIDTH;
    explicit BarcodeId(unsigned long value = 0) : value_(value) {}
    BarcodeId(unsigned long tile, unsigned long barcode, unsigned long cluster, unsigned long mismatches)
        : value_(((tile  & TILE_MASK) << TILE_SHIFT) |
                 ((barcode & BARCODE_MASK) << BARCODE_SHIFT) |
                 ((cluster & CLUSTER_MASK) << CLUSTER_SHIFT) |
                 ((mismatches & MISMATCHES_MASK) << MISMATCHES_SHIFT))
    {
        using boost::mpl::equal_to;
        using boost::mpl::int_;
        BOOST_MPL_ASSERT((equal_to<int_<57>, int_<MISMATCHES_WIDTH + CLUSTER_WIDTH + BARCODE_WIDTH + TILE_WIDTH> >));
        BOOST_MPL_ASSERT((equal_to<int_<57>, int_<TILE_WIDTH + TILE_SHIFT> >));
        using boost::format;
        using isaac::common::PreConditionException;
        if ((TILE_MASK < tile) |
            (BARCODE_MASK < barcode) |
            (CLUSTER_MASK < cluster) |
            (MISMATCHES_MASK < mismatches))
        {
            using boost::format;
            using isaac::common::PreConditionException;
            const format message = format(
                "BarcodeId(%ld, %ld, %ld, %ld): maximum values are (%ld, %ld, %ld, %ld)") %
                tile % barcode % cluster % mismatches %
                (unsigned long)TILE_MASK % (unsigned long)BARCODE_MASK % (unsigned long)CLUSTER_MASK % (unsigned long)MISMATCHES_MASK;
            BOOST_THROW_EXCEPTION(PreConditionException(message.str()));
        }
    }

    unsigned long getTile() const {return (value_ >> TILE_SHIFT) & TILE_MASK;}
    unsigned long getBarcode() const {return (value_ >> BARCODE_SHIFT) & BARCODE_MASK;}
    unsigned long getCluster() const {return (value_ >> CLUSTER_SHIFT) & CLUSTER_MASK;}
    unsigned long getMismatches() const {return (value_ >> MISMATCHES_SHIFT) & MISMATCHES_MASK;}
    unsigned long getTileBarcode() const {return (value_ >> BARCODE_SHIFT);}
    unsigned long getTileBarcodeCluster() const {return (value_ >> CLUSTER_SHIFT);}
    operator unsigned long() const {return value_;}
private:
    unsigned long value_;
};

inline std::ostream &operator<<(std::ostream &os, const BarcodeId &b)
{
    return os << "BarcodeId(" << b.getTile() << ":" << b.getBarcode() << ":" << b.getCluster() << ":" << b.getMismatches() << ")";
}

/**
 ** \brief Barcode bases with the information about the source cluster barcode mapping.
 **/
class Barcode
{
public:
    Barcode() : sequence_(0), barcodeId_(0) {}
    Barcode(Kmer sequence, BarcodeId barcodeId) : sequence_(sequence), barcodeId_(barcodeId) {}

    static Barcode constructFromTileBarcodeCluster(unsigned long tile, unsigned long barcode, unsigned long cluster)
    {
        return Barcode(0, BarcodeId(tile, barcode, cluster, 0));
    }

/*
    static Barcode constructFromSequenceBarcode(Kmer sequence, unsigned long barcode)
    {
        return Barcode(sequence, BarcodeId(0, barcode, 0));
    }
*/

    Barcode &operator=(const Barcode &barcode)
    {
        if (this != &barcode)
        {
            sequence_ = barcode.sequence_;
            barcodeId_ = barcode.barcodeId_;
        }
        return *this;
    }
    Kmer getSequence() const {return sequence_;}
    BarcodeId getBarcodeId() const {return barcodeId_;}
    unsigned long getTile() const {return barcodeId_.getTile();}
    unsigned long getBarcode() const {return barcodeId_.getBarcode();}
    unsigned long getCluster() const {return barcodeId_.getCluster();}
    unsigned long getMismatches() const {return barcodeId_.getMismatches();}
    void setSequence(Kmer bases) {sequence_ = bases;}
    void setBarcodeId(BarcodeId barcodeId) {barcodeId_ = barcodeId;}
private:
    Kmer sequence_;
    BarcodeId barcodeId_;
};

inline std::ostream &operator<<(std::ostream &os, const Barcode &barcode)
{
    return os << "Barcode(0x" << std::setfill('0') << std::setw(16) << std::hex << barcode.getSequence() <<
        "(" << bases(barcode.getSequence(), MAX_BARCODE_LENGTH) << ")," << barcode.getBarcodeId() << ")";
}

} // namespace demultiplexing
} // namespace isaac

#endif // #ifndef iSAAC_DEMULTIPLEXING_BARCODE_HH
