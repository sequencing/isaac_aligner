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
 ** \file SeedId.hh
 **
 ** \brief Identification of a seed encoding its tile, cluster, position and
 ** orientation.
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_SEED_ID_HH
#define iSAAC_ALIGNMENT_SEED_ID_HH

#include <cassert>
#include <iostream>
#include <boost/format.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/int.hpp>

#include "common/Exceptions.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Structured unique identifier of a seed.
 **
 ** From LSB to MSB:
 **   - reverse:  1 (            2)
 **   - seed   :  8 (          256)
 **   - cluster: 31 (2,147,483,648)
 **   - barcode: 12 (        4,096)
 **   - tile   : 12 (        4,096)
 **
 **
 ** The order of the fields is important as it will definet the natural order of
 ** the seed ids when sorted.
 **
 ** \note The implementation choice at the moment is to use a tile number that
 ** is the sequential index (0-based) of the tile in the input data set. This
 ** solution has been selected instead of keeping the actual tile and lane, to
 ** make easier (more flexible) to retrieve the metadata and components
 ** associated to each tile throughout the application.
 **
 ** \sa isaac::demultiplexing::BarcodeId
 **/
class SeedId
{
public:
    // width in bits for each field
    static const unsigned REVERSE_WIDTH = 1;
    static const unsigned SEED_WIDTH = 8;
    static const unsigned CLUSTER_WIDTH = 31;
    static const unsigned BARCODE_WIDTH = 12;
    static const unsigned TILE_WIDTH = 12;
    // masks for the values in each field
    static const unsigned long REVERSE_MASK = ~(~0UL<<REVERSE_WIDTH);
    static const unsigned long SEED_MASK = ~(~0UL<<SEED_WIDTH);
    static const unsigned long CLUSTER_MASK = ~(~0UL<<CLUSTER_WIDTH);
    static const unsigned long BARCODE_MASK = ~(~0UL<<BARCODE_WIDTH);
    static const unsigned long TILE_MASK = ~(~0UL<<TILE_WIDTH);
    // shifts in bits for each field
    static const unsigned REVERSE_SHIFT = 0;
    static const unsigned SEED_SHIFT = REVERSE_SHIFT + REVERSE_WIDTH;
    static const unsigned CLUSTER_SHIFT = SEED_SHIFT + SEED_WIDTH;
    static const unsigned BARCODE_SHIFT = CLUSTER_SHIFT + CLUSTER_WIDTH;
    static const unsigned TILE_SHIFT = BARCODE_SHIFT + BARCODE_WIDTH;
    explicit SeedId(unsigned long value = 0) : value_(value) {}
    SeedId(unsigned long tile, unsigned long barcode, unsigned long cluster, unsigned long seed, unsigned long reverse)
        : value_(((tile  & TILE_MASK) << TILE_SHIFT) |
                 ((barcode & BARCODE_MASK) << BARCODE_SHIFT) |
                 ((cluster & CLUSTER_MASK) << CLUSTER_SHIFT) |
                 ((seed & SEED_MASK) << SEED_SHIFT) |
                 ((reverse & REVERSE_MASK) << REVERSE_SHIFT))
    {
        using boost::mpl::equal_to;
        using boost::mpl::int_;
        BOOST_MPL_ASSERT((equal_to<int_<64>, int_<REVERSE_WIDTH + SEED_WIDTH + CLUSTER_WIDTH + BARCODE_WIDTH + TILE_WIDTH> >));
        BOOST_MPL_ASSERT((equal_to<int_<64>, int_<TILE_WIDTH + TILE_SHIFT> >));
        using boost::format;
        using isaac::common::PreConditionException;
        if ((TILE_MASK < tile) |
            (BARCODE_MASK < barcode) |
            (CLUSTER_MASK < cluster) |
            (SEED_MASK < seed) |
            (REVERSE_MASK < reverse))
        {
            using boost::format;
            using isaac::common::PreConditionException;
            const format message = format(
                "SeqId(%ld, %ld, %ld, %ld, %ld): maximum values are (%ld, %ld, %ld, %ld, %ld)") %
                tile % barcode % cluster % seed % reverse %
                (unsigned long)TILE_MASK % (unsigned long)BARCODE_MASK % (unsigned long)CLUSTER_MASK %
                (unsigned long)SEED_MASK % (unsigned long)REVERSE_MASK;
            BOOST_THROW_EXCEPTION(PreConditionException(message.str()));
        }
    }
    unsigned long getTile() const {return (value_ >> TILE_SHIFT) & TILE_MASK;}
    unsigned long getBarcode() const {return (value_ >> BARCODE_SHIFT) & BARCODE_MASK;}
    unsigned long getCluster() const {return (value_ >> CLUSTER_SHIFT) & CLUSTER_MASK;}
    unsigned long getTileBarcode() const {return (value_ >> BARCODE_SHIFT);}
    unsigned long getTileBarcodeCluster() const {return (value_ >> CLUSTER_SHIFT);}
    unsigned long getSeed() const {return (value_ >> SEED_SHIFT) & SEED_MASK;}
    bool isNSeedId() const {return SEED_MASK == getSeed();}
    void setNSeedId(const bool lowestSeed) {
        value_ &= ~(REVERSE_MASK << REVERSE_SHIFT);
        value_ |= (SEED_MASK << SEED_SHIFT) | (!lowestSeed << REVERSE_SHIFT);}
    bool isLowestNSeedId() const {return isNSeedId() && !isReverse();}
    unsigned long getReverse() const {return (value_ >> REVERSE_SHIFT) & REVERSE_MASK;}
    bool isReverse() const {return getReverse();}
    operator unsigned long() const {return value_;}
private:
    unsigned long value_;
};

/**
 * \brief A special SeedId that is greater than any normal seed id.
 */
const SeedId SMALLEST_N_SEED_ID(0, 0, 0, SeedId::SEED_MASK, false);

inline std::ostream &operator<<(std::ostream &os, const SeedId &s)
{
    return os << "SeedId(" << s.getTile() << ":" << s.getBarcode() << ":"<< s.getCluster() << ":" << s.getSeed() << ":" << s.getReverse() << ")";
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SEED_ID_HH
