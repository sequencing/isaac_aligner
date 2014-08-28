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
 ** \file Kmer.hh
 **
 ** General definitions and tools for handling k-mers.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_OLIGO_K_MER_HH
#define iSAAC_OLIGO_K_MER_HH

#include <string>

#include "oligo/Nucleotides.hh"

namespace isaac
{
namespace oligo
{

typedef __uint128_t LongKmerType;
BOOST_STATIC_ASSERT(16 == sizeof(LongKmerType));

typedef unsigned long KmerType;
BOOST_STATIC_ASSERT(8 == sizeof(KmerType));

typedef unsigned ShortKmerType;
BOOST_STATIC_ASSERT(4 == sizeof(ShortKmerType));


static const unsigned BITS_PER_BASE = 2;
static const unsigned BITS_PER_BASE_MASK = 3;

template <typename KmerT>
struct KmerTraits
{
    static const unsigned KMER_BASES = sizeof(KmerT) * 8 / BITS_PER_BASE;
    static const unsigned KMER_BITS = BITS_PER_BASE * KMER_BASES;
};

template <typename KmerT>
std::string bases(KmerT kmer)
{
    return bases<BITS_PER_BASE>(kmer, KmerTraits<KmerT>::KMER_BASES);
}

template <typename KmerT>
inline std::string reverseBases(KmerT kmer)
{
    std::string s = bases(~kmer);
    std::reverse(s.begin(), s.end());
    return s;
}

template <typename KmerT>
inline std::ostream & operator <<(std::ostream &os, const oligo::Bases<BITS_PER_BASE, KmerT> &bases)
{
    return printBases(os, bases);
}

template <typename KmerT>
inline std::ostream & operator <<(std::ostream &os, const oligo::ReverseBases<BITS_PER_BASE, KmerT> &bases)
{
    return printReverseBases(os, bases);
}

template <typename KmerT>
struct TraceKmer
{
    KmerT kmer_;
    TraceKmer(KmerT kmer) : kmer_(kmer) {}
};

template <typename KmerT>
TraceKmer<KmerT> traceKmer(const KmerT &kmer)
{
    return TraceKmer<KmerT>(kmer);
}

template <typename T> std::ostream &traceKmerValue(
	std::ostream &os, const T &kmer)
{
    return os << "0x" << std::hex << kmer;
}

template <> inline std::ostream &traceKmerValue<__uint128_t>(
	std::ostream &os, const __uint128_t &kmer)
{
    const unsigned long lobits(kmer);
    const unsigned long hibits(kmer >> 64);
    if (hibits)
    {
        return os << "0x" << std::hex << hibits << lobits;
    }
    return os << "0x" << std::hex << lobits;
}

/**
 * \brief special type to ensure that serialization of Kmer compiles regardless of 
 *        how many bits it is defined for. Currently std::osream << __uint128_t does not
 *        compile well
 */
template <typename KmerT>
inline std::ostream & operator << (std::ostream &os, const TraceKmer<KmerT> &kmer)
{
    return traceKmerValue(os, kmer.kmer_);
}

} // namespace oligo
} // namespace isaac

#endif // #ifndef iSAAC_OLIGO_K_MER_HH
