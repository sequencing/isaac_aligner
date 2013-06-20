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
 ** \file Permutate.cpp
 **
 ** \brief See Permutate.hh.
 **
 ** \author Come Raczy
 **/

#include <iostream>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "oligo/Permutate.hh"

namespace isaac
{
namespace oligo
{

Permutate::Permutate(const unsigned blockLength, const std::vector<unsigned> &from, const std::vector<unsigned> &to)
    : blockLength_(blockLength)
    , count_(from.size())
    , order_(encode(from, to))
    , absoluteReverseOrder_(encode(to))\
    , from_(from)
    , to_(to)
{
}

std::string Permutate::toString() const
{
    static std::string S = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    std::string ret = "from ";
    BOOST_FOREACH(unsigned v, from_)
    {
        ret.push_back(v < S.size() ? S[v] : '?');
    }
    ret += " to ";
    BOOST_FOREACH(unsigned v, to_)
    {
        ret.push_back(v < S.size() ? S[v] : '?');
    }
    return ret;
}

unsigned long Permutate::encode(const std::vector<unsigned> &from, const std::vector<unsigned> &to)
{
    //std::cerr << "\nENCODE" << std::endl;
    assert(from.size() == to.size());
    assert(from.size() <= 16);
    unsigned long ret = 0;
    std::vector<bool> checkFrom(from.size(), false);
    std::vector<bool> checkTo(to.size(), false);
    for (unsigned origin = 0; from.size() > origin; ++origin)
    {
        assert(!checkFrom[origin]);
        assert(!checkTo[origin]);
        checkFrom[origin] = true;
        checkTo[origin] = true;
        const std::vector<unsigned>::const_iterator found = std::find(to.begin(), to.end(), from[origin]);
        assert(to.end() != found);
        const unsigned target = found - to.begin();
        assert(target < from.size());
        ret <<= ENCODING_BITS;
        ret |= target;
        //std::cerr << (boost::format("%d:%d:%d:%08x") % origin % target % from[origin] % ret).str() << std::endl;
    }

    assert(0 == std::count(checkFrom.begin(), checkFrom.end(), false));
    assert(0 == std::count(checkTo.begin(), checkTo.end(), false));
    return ret;
}

unsigned long Permutate::encode(const std::vector<unsigned> &to)
{
    unsigned i = 0;
    using boost::lambda::var;
    using namespace boost::assign_detail;
    const std::vector<unsigned> origin = generic_list<unsigned>().repeat_fun(to.size(), var(i)++);
    return encode(to, origin);
}

void buildPermutationList(
    const std::vector<unsigned> &prefix,
    const std::vector<unsigned> &suffix,
    const unsigned n,
    std::vector<std::vector<unsigned> > &permutationList)
{
    assert(2 * n == prefix.size() + suffix.size());
    if (prefix.size() == n)
    {
        permutationList.push_back(prefix);
        std::vector<unsigned> &back = permutationList.back();
        back.insert(back.end(), suffix.begin(), suffix.end());
    }
    else
    {
        for (unsigned i = 0; suffix.size() > i; ++i)
        {
            if (prefix.empty() || suffix[i] > prefix.back())
            {
                std::vector<unsigned> newPrefix = prefix;
                newPrefix.push_back(suffix[i]);
                std::vector<unsigned> newSuffix = suffix;
                newSuffix.erase(newSuffix.begin() + i);
                buildPermutationList(newPrefix, newSuffix, n, permutationList);
            }
        }
    }
}

template <typename KmerT>
std::vector<Permutate> getPermutateList(const unsigned errorCount)
{
    const unsigned blocksCount = 2 * errorCount;
    const unsigned blockLength = oligo::KmerTraits<KmerT>::KMER_BASES / blocksCount;
    assert(blocksCount * blockLength == oligo::KmerTraits<KmerT>::KMER_BASES);
    const std::vector<unsigned> prefix;
    using boost::assign_detail::generic_list;
    using boost::lambda::var;
    unsigned i = 0;
    const std::vector<unsigned> suffix = generic_list<unsigned>().repeat_fun(2 * errorCount, var(i)++);
    std::vector<std::vector<unsigned> > permutationList;
    buildPermutationList(prefix, suffix, errorCount, permutationList);
    std::vector<oligo::Permutate> permutateList;
    std::vector<std::vector<unsigned> >::const_iterator from = permutationList.begin();
    for (std::vector<std::vector<unsigned> >::const_iterator to = from; permutationList.end() != to; ++to)
    {
        permutateList.push_back(Permutate(blockLength, *from, *to));
        from = to;
    }
    return permutateList;
}


template<typename KmerT>
KmerT Permutate::transform(const KmerT kmer, const unsigned long order) const
{
    ISAAC_ASSERT_MSG(blockLength_ * count_ <= oligo::KmerTraits<KmerT>::KMER_BASES, "Permutation is incompatible with Kmer type");
    ISAAC_ASSERT_MSG(count_ <= (1 << ENCODING_BITS), "Block count should fit in " << ENCODING_BITS << " bits");

    const unsigned blockBits = BITS_PER_BASE * blockLength_;
    const KmerT blockMask = ~((~KmerT(0)) << blockBits);
    //std::cerr << (boost::format("\n%d:%016x:%016x:%d:%016x") % blockBits % blockMask % order_ % count_ %kmer).str() << std::endl;
    KmerT ret = 0;
    for (unsigned origin = 0; count_ > origin; ++origin)
    {
        const unsigned orderShift = (count_ - origin - 1U) * ENCODING_BITS;
        const unsigned target = (order >> orderShift) & ENCODING_MASK;
        const unsigned kmerOriginShift = (count_ - origin - 1U) * blockBits;
        const unsigned kmerTargetShift = (count_ - target - 1U) * blockBits;
        ret |= (((kmer >> kmerOriginShift) & blockMask) << kmerTargetShift);
        //std::cerr << (boost::format("    %016x:%d:%d:%d:%d") % ret % origin % target % kmerOriginShift % kmerTargetShift).str() << std::endl;
    }
    return ret;
}

template oligo::ShortKmerType Permutate::transform<oligo::ShortKmerType>(const oligo::ShortKmerType kmer, const unsigned long order) const;
template oligo::KmerType Permutate::transform<oligo::KmerType>(const oligo::KmerType kmer, const unsigned long order) const;
template oligo::LongKmerType Permutate::transform<oligo::LongKmerType>(const oligo::LongKmerType kmer, const unsigned long order) const;

template std::vector<Permutate> getPermutateList<oligo::ShortKmerType>(const unsigned errorCount);
template std::vector<Permutate> getPermutateList<oligo::KmerType>(const unsigned errorCount);
template std::vector<Permutate> getPermutateList<oligo::LongKmerType>(const unsigned errorCount);

} //namespace oligo
} // namespace isaac
