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
 ** \file Seed.hh
 **
 ** \brief Definition od a seed by its k-mer and identifier (SeedId).
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_SEED_HH
#define iSAAC_ALIGNMENT_SEED_HH

#include <cassert>
#include <iostream>
#include <boost/format.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/int.hpp>

#include "oligo/Kmer.hh"
#include "alignment/SeedId.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Structured unique identifier of a seed.
 **/
template <typename KmerT>
class Seed
{
public:
    Seed() : kmer_(0), seedId_(0) {}
    Seed(KmerT kmer, SeedId seedId) : kmer_(kmer), seedId_(seedId) {}
    Seed(const Seed &seed) : kmer_(seed.kmer_), seedId_(seed.seedId_) {}
    Seed &operator=(const Seed &seed)
    {
        if (this != &seed)
        {
            kmer_ = seed.kmer_;
            seedId_ = seed.seedId_;
        }
        return *this;
    }
    KmerT &kmer() {return kmer_;}
    KmerT getKmer() const {return kmer_;}
    SeedId getSeedId() const {return seedId_;}
    unsigned long getTile() const {return seedId_.getTile();}
    unsigned long getBarcode() const {return seedId_.getBarcode();}
    unsigned long getCluster() const {return seedId_.getCluster();}
    unsigned long getSeedIndex() const {return seedId_.getSeed();}
    bool isNSeed() const {return seedId_.isNSeedId();}
    bool isLowestNSeed() const {return seedId_.isLowestNSeedId();}
    /**
     * \brief produces an N-seed
     *        Not all N-seeds are equal. The ones that have been built out of
     *        seed with index 0, have their reverse bit set to false. This allows distinct behavior when
     *        storing no-matches in MatchFinder
     */
    void makeNSeed(bool lowestNSeed) {kmer_ = ~KmerT(0); seedId_.setNSeedId(lowestNSeed);}
    bool isReverse() const {return getSeedId().isReverse();}
    void setKmer(KmerT kmer) {kmer_ = kmer;}
    void setSeedId(SeedId seedId) {seedId_ = seedId;}
private:
    KmerT kmer_;
    SeedId seedId_;
};

template <typename KmerT>
inline Seed<KmerT> makeNSeed(unsigned long tile, unsigned long barcode, unsigned long cluster, bool lowestSeedId)
{
    return Seed<KmerT>(~KmerT(0), SeedId(tile, barcode, cluster, SeedId::SEED_MASK, !lowestSeedId));
}

template <typename KmerT>
inline bool orderByKmerSeedIndex(const Seed<KmerT> &lhs, const Seed<KmerT> &rhs)
{
    // IMPORTANT!!!: The match finder relies on N-seeds to be at the end of the seed list after sorting by kmer.
    // N-seeds are assigned the highest possible seed index value by seed loader
    return (lhs.getKmer() < rhs.getKmer()) || (lhs.getKmer() == rhs.getKmer() && lhs.getSeedId().getSeed() < rhs.getSeedId().getSeed());
}

template <typename KmerT>
inline std::ostream &operator<<(std::ostream &os, const Seed<KmerT> &seed)
{
    return os << "Seed(" << oligo::Bases<oligo::BITS_PER_BASE, KmerT>(seed.getKmer(), oligo::KmerTraits<KmerT>::KMER_BASES) <<
        "(" << oligo::ReverseBases<oligo::BITS_PER_BASE, KmerT>(seed.getKmer(), oligo::KmerTraits<KmerT>::KMER_BASES) << ")" <<
        "," << seed.getSeedId() << ")";
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SEED_HH
