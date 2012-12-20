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
 ** \file Nucleotides.hh
 **
 ** General tools and definitions to manipulate nucleotides.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_OLIGO_NUCLEOTIDES_HH
#define iSAAC_OLIGO_NUCLEOTIDES_HH

#include <vector>
#include <boost/array.hpp>
#include <boost/assign.hpp>

namespace isaac
{
namespace oligo
{

// valid ones are 0(A) 1(C) 2(G) and 3(T).
// for data, invalidOligo represents N
// for reference, invalidOligo indicates any non-ACGT base value
static const unsigned int invalidOligo = 4;

inline std::vector<unsigned int> getTranslator(const bool withN = false, const unsigned defaultValue = invalidOligo)
{   
    std::vector<unsigned int> translator(256, defaultValue);
    translator['a'] = 0;
    translator['A'] = 0;
    translator['c'] = 1;
    translator['C'] = 1;
    translator['g'] = 2;
    translator['G'] = 2;
    translator['t'] = 3;
    translator['T'] = 3;
    if (withN)
    {
        translator['n'] = invalidOligo;
        translator['N'] = invalidOligo;
    }
    return translator;
}

extern const std::vector<unsigned int> defaultTranslator;

inline unsigned int getValue(const char base)
{   
    static const std::vector<unsigned int> translator = getTranslator();
    assert(0 < base);
    assert(base < static_cast<long>(translator.size()));
    return translator[base];
}

static const boost::array<char, 5> allBases = boost::assign::list_of('A')('C')('G')('T')('N');
inline char getBase(unsigned int base, const bool upperCase = true)
{
    base = std::min(base, static_cast<unsigned int>(allBases.size()));    
    if (upperCase)
    {
        return allBases[base];
    }
    else
    {
        return allBases[base] + ('a' - 'A');
    }
}

inline char getUppercaseBase(unsigned int base)
{
    return getBase(base, true);
}

inline char getUppercaseBaseFromBcl(unsigned char bcl)
{
    return getBase(bcl & 0x03, true);
}

static const boost::array<char, 5> allReverseBases = boost::assign::list_of('T')('G')('C')('A')('N');
inline char getReverseBase(unsigned int base, const bool upperCase = true)
{
//    return base > 3 ? getBase(base, upperCase) : getBase((~base) & 3, upperCase);
    base = std::min(base, static_cast<unsigned int>(allReverseBases.size()));
    if (upperCase)
    {
        return allReverseBases[base];
    }
    else
    {
        return allReverseBases[base] + ('a' - 'A');
    }
}

inline char getReverseBase(const char base)
{
    switch (base)
    {
    case 'a': return 't';
    case 'A': return 'T';
    case 'c': return 'g';
    case 'C': return 'G';
    case 'g': return 'c';
    case 'G': return 'C';
    case 't': return 'a';
    case 'T': return 'A';
    case 'n': return 'n';
    default : return 'N';
    }
}

inline bool isBclN(const char bclByte)
{
    return !(bclByte & 0xfc);
}

/**
 * \brief reverse-complements the base bits of a bcl byte
 *
 * \returns 0 for 0, reverse-complemented lower bits with quality bits unchanged
 */
inline char getReverseBcl(const char bcl)
{
    return !isBclN(bcl) ? (bcl & 0xfc) | (0x3 - (bcl & 0x03)) : 0;
}

template <int bitsPerBase, typename KmerT>
std::string bases(KmerT kmer, unsigned kmerLength)
{
    static const KmerT kmerMask = ~(KmerT(0)) >> (sizeof(KmerT) * 8 - bitsPerBase);

    std::string s;
    while (kmerLength--)
    {
        s.push_back(isaac::oligo::getBase((kmer >> bitsPerBase * kmerLength) & kmerMask));
    }
    return s;
}

template <unsigned bitsPerBase, typename KmerT>
struct Bases
{
    const KmerT kmer_;
    const unsigned kmerLength_;
    static const unsigned bitsPerBase_ = bitsPerBase;
    static const KmerT kmerMask_ = ~(KmerT(0)) >> (sizeof(KmerT) * 8 - bitsPerBase);
    Bases(const KmerT kmer, unsigned kmerLength) : kmer_(kmer), kmerLength_(kmerLength)
    {

    }
};

template <unsigned bitsPerBase, typename KmerT>
struct ReverseBases
{
    const KmerT kmer_;
    const unsigned kmerLength_;
    static const unsigned bitsPerBase_ = bitsPerBase;
    static const KmerT kmerMask_ = ~(KmerT(0)) >> (sizeof(KmerT) * 8 - bitsPerBase);
    ReverseBases(const KmerT kmer, unsigned kmerLength) : kmer_(kmer), kmerLength_(kmerLength)
    {

    }
};

template <typename BasesT>
std::ostream & printBases(std::ostream &os, const BasesT &b)
{
    unsigned pos = b.kmerLength_;
    while (pos--)
    {
        os << isaac::oligo::getBase((b.kmer_ >> BasesT::bitsPerBase_ * pos) & BasesT::kmerMask_);
    }
    return os;
}

template <typename BasesT>
std::ostream & printReverseBases(std::ostream &os, const BasesT &b)
{
    for (unsigned pos = 0; b.kmerLength_ > pos; ++pos)
    {
        os << isaac::oligo::getReverseBase(unsigned((b.kmer_ >> BasesT::bitsPerBase_ * pos) & BasesT::kmerMask_));
    }
    return os;
}

inline std::string bclToString(const unsigned char *basesIterator, unsigned length)
{
    std::string ret;
    while(length--)
    {
        ret += oligo::getBase(0x3 & *basesIterator++);
    }
    return ret;
}

template <typename FwdIteratorT>
unsigned long pack32BclBases(FwdIteratorT bcl)
{
    unsigned long ret = 0;
    ret |= (0x3ul & *bcl++);
    ret |= (0x3ul & *bcl++) <<  2;
    ret |= (0x3ul & *bcl++) <<  4;
    ret |= (0x3ul & *bcl++) <<  6;
    ret |= (0x3ul & *bcl++) <<  8;
    ret |= (0x3ul & *bcl++) << 10;
    ret |= (0x3ul & *bcl++) << 12;
    ret |= (0x3ul & *bcl++) << 14;
    ret |= (0x3ul & *bcl++) << 16;
    ret |= (0x3ul & *bcl++) << 18;
    ret |= (0x3ul & *bcl++) << 20;
    ret |= (0x3ul & *bcl++) << 22;
    ret |= (0x3ul & *bcl++) << 24;
    ret |= (0x3ul & *bcl++) << 26;
    ret |= (0x3ul & *bcl++) << 28;
    ret |= (0x3ul & *bcl++) << 30;
    ret |= (0x3ul & *bcl++) << 32;
    ret |= (0x3ul & *bcl++) << 34;
    ret |= (0x3ul & *bcl++) << 36;
    ret |= (0x3ul & *bcl++) << 38;
    ret |= (0x3ul & *bcl++) << 40;
    ret |= (0x3ul & *bcl++) << 42;
    ret |= (0x3ul & *bcl++) << 44;
    ret |= (0x3ul & *bcl++) << 46;
    ret |= (0x3ul & *bcl++) << 48;
    ret |= (0x3ul & *bcl++) << 50;
    ret |= (0x3ul & *bcl++) << 52;
    ret |= (0x3ul & *bcl++) << 54;
    ret |= (0x3ul & *bcl++) << 56;
    ret |= (0x3ul & *bcl++) << 58;
    ret |= (0x3ul & *bcl++) << 60;
    ret |= (0x3ul & *bcl)   << 62;
    return ret;
}

} // namespace oligo
} // namespace isaac

#endif // #ifndef iSAAC_OLIGO_NUCLEOTIDES_HH
