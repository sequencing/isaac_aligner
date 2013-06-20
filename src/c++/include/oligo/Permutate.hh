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
 ** \file Permutate.hh
 **
 ** \brief Utility class to permutate blocks in a kmer.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_OLIGO_PERMUTATE_HH
#define iSAAC_OLIGO_PERMUTATE_HH

#include <vector>
#include <string>

#include <boost/lambda/lambda.hpp>

#include "oligo/Kmer.hh"

namespace isaac
{
namespace oligo
{

/**
 ** \brief Utility class to reorder the blocks in a kmer, from an origin permutation
 ** to a target permutation.
 **
 ** The kmer is split into 'count' (<=16) blocks, each of them at positions [0,
 ** count). The number of bases in each block is constant.
 **/
class Permutate
{
public:
    /// Create a permutation from a list of positions. The un-permuted order is [0,count),
    Permutate(unsigned blockLength, const std::vector<unsigned> &from, const std::vector<unsigned> &to);

    /// Apply the permutation to the kmer
    template <typename KmerT>
    KmerT operator()(const KmerT kmer) const
    {
        return transform(kmer, order_);
    }

    /// Reorder tyhe blocks in their natural sequence (0, 1, 2...)
    template <typename KmerT>
    KmerT reorder(const KmerT kmer) const
    {
        return transform(kmer, absoluteReverseOrder_);
    }
    /// Shows a readable representation of the permutation
    std::string toString() const;
private:
   /// the length of the blocks
    unsigned blockLength_;
    /// the count of blocks
    unsigned count_;
    /// the encoded order of the blocks relatively to 'from'
    unsigned long order_;
    /// the encoded absolute order of the blocks relatively to the 'to' order
    unsigned long absoluteReverseOrder_;
    ///
    std::vector<unsigned> from_;
    std::vector<unsigned> to_;
    /**
     ** \brief encode the re-ordering of the blocks from a permutation to an other.
     **
     ** Example
     **
     ** Assuming 4 blocks ABCD, numbered respectively 0, 1, 2 and 3. If the
     ** original permutation is ABDC and the targeted permutation ACBD, the
     ** 'from' vector would be [0, 1, 3, 2] and the 'to' vector [0, 2, 1,
     ** 3]. The encoded value would be 0x0231 (the block at position 0 stays at
     ** position 0, position 1 goes to position 2, position 2 goes to position 3
     ** and position 3 goes to position 1).
     **
     ** \param from the list of block number in the original permutation
     **
     ** \param to the list of block numbers in the targetted permutation
     **
     ** \return the encoded transformation from/to
     **/
    unsigned long encode(const std::vector<unsigned> &from, const std::vector<unsigned> &to);
    /// encode the re-ordering of the blocks from the natural order (0, 1, 2, ...)
    unsigned long encode(const std::vector<unsigned> &to);
    /// apply the permutation encoded in the given order to the kmer
    template<typename KmerT>
    KmerT transform(const KmerT kmer, const unsigned long order) const;

    static const unsigned ENCODING_BITS = 4;
    static const unsigned long ENCODING_MASK = 0x0F;
};

/**
 ** \brief Generate the list of permutations for a given number of errors
 **
 ** \param errorCount the number of errors to support
 **
 ** \return the list of Permutate components in the order where they should be
 ** applied (starting from the natural order).
 **/
template <typename KmerT>
std::vector<Permutate> getPermutateList(const unsigned errorCount);


} // namespace oligo
} // namespace isaac

#endif // #ifndef iSAAC_OLIGO_PERMUTATE_HH
