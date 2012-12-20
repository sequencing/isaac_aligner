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

typedef unsigned long Kmer;
typedef unsigned int  NMask;
//typedef std::pair<Kmer, reference::ReferencePosition> Pair;
const unsigned BITS_PER_BASE = 2;
const unsigned int kmerLength = 32;
const unsigned int kmerBitLength = BITS_PER_BASE * kmerLength;

inline std::string bases(Kmer kmer)
{
    return bases<BITS_PER_BASE>(kmer, kmerLength);
}

inline std::string reverseBases(Kmer kmer)
{
    std::string s = bases(~kmer);
    std::reverse(s.begin(), s.end());
    return s;
}

inline std::ostream & operator <<(std::ostream &os, const oligo::Bases<BITS_PER_BASE, Kmer> &bases)
{
    return printBases(os, bases);
}

inline std::ostream & operator <<(std::ostream &os, const oligo::ReverseBases<BITS_PER_BASE, Kmer> &bases)
{
    return printReverseBases(os, bases);
}

} // namespace oligo
} // namespace isaac

#endif // #ifndef iSAAC_OLIGO_K_MER_HH
