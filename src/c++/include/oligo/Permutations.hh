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
 ** \file Permutations.hh
 **
 ** Permutation operations on oligomers.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_OLIGO_PERMUTATIONS_HH
#define iSAAC_OLIGO_PERMUTATIONS_HH

#include <iostream>
#include <utility>
#include <string>
#include <map>

#include <boost/foreach.hpp>
#include <boost/assign.hpp>

#include "common/Exceptions.hh"
#include "oligo/Kmer.hh"
#include "oligo/Mask.hh"
#include "reference/ReferenceKmer.hh"
#include "common/SystemCompatibility.hh"

namespace isaac
{
namespace oligo
{
const unsigned int blockLength = kmerLength / 4;
const unsigned int blockBitLength = 2 * blockLength;
//const Kmer blockMask = 0xFFFFUL;

typedef Kmer(*Permutation)(Kmer);
typedef std::pair<Permutation, std::string> NamedPermutation;

inline Kmer pABCD2ABCD(const Kmer kmer);

inline Kmer pBCDA2ABCD(const Kmer kmer);
inline Kmer pCDAB2ABCD(const Kmer kmer);
inline Kmer pACBD2ABCD(const Kmer kmer);
inline Kmer pBDAC2ABCD(const Kmer kmer);
inline Kmer pADBC2ABCD(const Kmer kmer);

inline Kmer pABCD2BCDA(const Kmer kmer);
inline Kmer pABCD2CDAB(const Kmer kmer);
inline Kmer pABCD2ACBD(const Kmer kmer);
inline Kmer pABCD2BDAC(const Kmer kmer);
inline Kmer pABCD2ADBC(const Kmer kmer);

inline Kmer pBCDA(const Kmer kmer);
inline Kmer pCDAB(const Kmer kmer);
inline Kmer pACBD(const Kmer kmer);
inline Kmer pBDAC(const Kmer kmer);
inline Kmer pADBC(const Kmer kmer);

const std::vector<NamedPermutation> permutations = boost::assign::list_of
    (NamedPermutation(&pBCDA, "BCDA"))
    (NamedPermutation(&pCDAB, "CDAB"))
    (NamedPermutation(&pACBD, "ACBD"))
    (NamedPermutation(&pBDAC, "BDAC"))
    (NamedPermutation(&pADBC, "ADBC"))
    ;

const std::vector<std::string> permutationNameList = boost::assign::list_of
    ("ABCD")("ACBD")("ADBC")("BCDA")("BDAC")("CDAB");

// returns direct revese permutation from the specified name to ABCD
inline Permutation getReversePermutation(const std::string &sourcePermutation)
{
    const std::map<std::string, Permutation > directPermutations = boost::assign::map_list_of
        ("ABCD", pABCD2ABCD)
        ("BCDA", pBCDA2ABCD)
        ("CDAB", pCDAB2ABCD)
        ("ACBD", pACBD2ABCD)
        ("BDAC", pBDAC2ABCD)
        ("ADBC", pADBC2ABCD)
        ;

    const std::map<std::string, Permutation >::const_iterator cit(directPermutations.find(sourcePermutation));
    if (directPermutations.end() == cit)
    {
        BOOST_THROW_EXCEPTION(common::InvalidParameterException("Unknown permutation name: " + sourcePermutation));
    }
    return cit->second;
}

// returns direct permutation from ABCD into the specified name
inline Permutation getPermutation(const std::string &targetPermutation)
{
    const std::map<std::string, Permutation > directPermutations = boost::assign::map_list_of
        ("ABCD", pABCD2ABCD)
        ("BCDA", pABCD2BCDA)
        ("CDAB", pABCD2CDAB)
        ("ACBD", pABCD2ACBD)
        ("BDAC", pABCD2BDAC)
        ("ADBC", pABCD2ADBC)
        ;

    const std::map<std::string, Permutation >::const_iterator cit(directPermutations.find(targetPermutation));
    if (directPermutations.end() == cit)
    {
        BOOST_THROW_EXCEPTION(common::InvalidParameterException("Unknown permutation name: " + targetPermutation));
    }
    return cit->second;
}

inline void permuteBlocks(Permutation permutation, std::vector<reference::ReferenceKmer> &reference)
{
    std::cerr << "Permuting " << reference.size() << " " << kmerLength << "-mers" << std::endl;
    const clock_t start = clock();
    BOOST_FOREACH(reference::ReferenceKmer &referenceKmer, reference)
    {
        referenceKmer.setKmer(permutation(referenceKmer.getKmer()));
    }
    std::cerr << "... done in " << (clock() - start) / 1000 << "ms" << std::endl;
}

// ABCD
inline Kmer pABCD2ABCD(const Kmer kmer)
{
    return kmer;
}

// BCDA
inline Kmer pBCDA(const Kmer kmer)
{
    return (kmer << blockBitLength) | (kmer >> (3 * blockBitLength));
}

inline Kmer pABCD2BCDA(const Kmer kmer)
{
    return (kmer << blockBitLength) | (kmer >> (3 * blockBitLength));
}

inline Kmer pBCDA2ABCD(const Kmer kmer)
{
    return (kmer >> blockBitLength) | (kmer << (3 * blockBitLength));
}

// CDAB
inline Kmer pCDAB(const Kmer kmer)
{
    return (kmer << blockBitLength) | (kmer >> (3 * blockBitLength));
}

inline Kmer pABCD2CDAB(const Kmer kmer)
{
    return pCDAB2ABCD(kmer);
}

inline Kmer pCDAB2ABCD(const Kmer kmer)
{
    return (kmer >> 2 * blockBitLength) | (kmer << (2 * blockBitLength));
}

// ACBD
inline Kmer pACBD(const Kmer kmer)
{
    const Kmer ACOO = ((kmer & 0xFFFF0000UL) | (kmer >> (3 * blockBitLength))) << (2 * blockBitLength);
    const Kmer OOBD = ((kmer << (3 * blockBitLength)) | (kmer & 0xFFFF0000FFFFUL)) >> (2 * blockBitLength);
    return ACOO | OOBD;
}

inline Kmer pABCD2ACBD(const Kmer kmer)
{
    return  pACBD2ABCD(kmer);
}

inline Kmer pACBD2ABCD(const Kmer kmer)
{
    return  (kmer                   & 0xFFFF00000000FFFFUL) |
            (kmer >> blockBitLength & 0x00000000FFFF0000UL) |
            (kmer << blockBitLength & 0x0000FFFF00000000UL);
}

// BDAC
inline Kmer pBDAC(const Kmer kmer)
{
    return (kmer << (2 * blockBitLength)) | (kmer >> (2 * blockBitLength));
}

inline Kmer pABCD2BDAC(const Kmer kmer)
{
    return  (kmer <<     blockBitLength & 0xFFFF000000000000UL) | //move B
            (kmer << 2 * blockBitLength & 0x0000FFFF00000000UL) | //move D
            (kmer >> 2 * blockBitLength & 0x00000000FFFF0000UL) | //move A
            (kmer >>     blockBitLength & 0x000000000000FFFFUL);  //move C
}

inline Kmer pBDAC2ABCD(const Kmer kmer)
{
    return  (kmer << 2 * blockBitLength & 0xFFFF000000000000UL) | //move A
            (kmer >>     blockBitLength & 0x0000FFFF00000000UL) | //move B
            (kmer <<     blockBitLength & 0x00000000FFFF0000UL) | //move C
            (kmer >> 2 * blockBitLength & 0x000000000000FFFFUL);  //move D
}

//ADBC
inline Kmer pADBC(const Kmer kmer)
{
    const Kmer ADOO = ((kmer & 0xFFFF0000UL) | ((kmer >> (2 * blockBitLength)) & 0xFFFFUL)) << (2 * blockBitLength);
    const Kmer OOBC = ((kmer >> (3 * blockBitLength)) << blockBitLength) | (kmer & 0xFFFFUL);
    return ADOO | OOBC;
}

inline Kmer pABCD2ADBC(const Kmer kmer)
{
    return  (kmer                       & 0xFFFF000000000000UL) | //keep A
            (kmer >>     blockBitLength & 0x00000000FFFFFFFFUL) | //move BC
            (kmer << 2 * blockBitLength & 0x0000FFFF00000000UL);  //move D
}

inline Kmer pADBC2ABCD(const Kmer kmer)
{
    return  (kmer                       & 0xFFFF000000000000UL) | //keep A
            (kmer <<     blockBitLength & 0x0000FFFFFFFF0000UL) | //move BC
            (kmer >> 2 * blockBitLength & 0x000000000000FFFFUL);  //move D
}

} // namespace oligo
} // namespace isaac

#endif // #ifndef iSAAC_OLIGO_PERMUTATIONS_HH

