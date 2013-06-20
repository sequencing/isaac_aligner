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
 **/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/math/tools/big_constant.hpp>
using namespace std;

#include "RegistryName.hh"
#include "testPermutate.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestPermutate, registryName("Permutate"));

void TestPermutate::setUp()
{
}

void TestPermutate::tearDown()
{
}

void TestPermutate::testFourBlocks()
{
    using namespace isaac;
    using boost::assign::list_of;
    const unsigned long kmer = 0xFEDCBA9876543210UL;
    const unsigned long abcd = 0xFEDCBA9876543210UL;
    const unsigned long adbc = 0xFEDC3210BA987654UL;
    const unsigned long dbca = 0x3210BA987654FEDCUL;
    const unsigned blockLength = 8;
    const std::vector<unsigned> ABCD = list_of(0)(1)(2)(3);
    const std::vector<unsigned> ADBC = list_of(0)(3)(1)(2);
    const std::vector<unsigned> DBCA = list_of(3)(1)(2)(0);
    const oligo::Permutate ABCD_ABCD(blockLength, ABCD, ABCD);
    const oligo::Permutate ABCD_ADBC(blockLength, ABCD, ADBC);
    const oligo::Permutate ABCD_DBCA(blockLength, ABCD, DBCA);
    const oligo::Permutate ADBC_DBCA(blockLength, ADBC, DBCA);
    CPPUNIT_ASSERT_EQUAL(ABCD_ABCD(kmer), kmer);
    CPPUNIT_ASSERT_EQUAL(ABCD_ADBC(kmer), 0xFEDC3210BA987654UL);
    CPPUNIT_ASSERT_EQUAL(ABCD_DBCA(kmer), 0x3210BA987654FEDCUL);
    CPPUNIT_ASSERT_EQUAL(ADBC_DBCA(kmer), 0xBA9876543210FEDCUL);
    CPPUNIT_ASSERT_EQUAL(ABCD_ABCD(abcd), abcd);
    CPPUNIT_ASSERT_EQUAL(ABCD_ADBC(abcd), adbc);
    CPPUNIT_ASSERT_EQUAL(ABCD_DBCA(abcd), dbca);
    CPPUNIT_ASSERT_EQUAL(ADBC_DBCA(adbc), dbca);
    CPPUNIT_ASSERT_EQUAL(ABCD_ABCD.reorder(abcd), abcd);
    CPPUNIT_ASSERT_EQUAL(ABCD_ADBC.reorder(adbc), abcd);
    CPPUNIT_ASSERT_EQUAL(ABCD_DBCA.reorder(dbca), abcd);
    CPPUNIT_ASSERT_EQUAL(ADBC_DBCA.reorder(dbca), abcd);
}

void TestPermutate::testEightBlocks()
{
    using namespace isaac;
    using boost::assign::list_of;
    // A  B  C  D  E  F  G  H
    //FE DC BA 98 76 54 32 10
    const unsigned long kmer = 0xFEDCBA9876543210UL;
    const unsigned long abcdefgh = 0xFEDCBA9876543210UL;
    const unsigned long adbcefgh = 0xFE98DCBA76543210UL;
    const unsigned long dghbcafe = 0x983210DCBAFE5476UL;
    const unsigned blockLength = 4;
    const std::vector<unsigned> ABCDEFGH = list_of(0)(1)(2)(3)(4)(5)(6)(7);
    const std::vector<unsigned> ADBCEFGH = list_of(0)(3)(1)(2)(4)(5)(6)(7);
    const std::vector<unsigned> DGHBCAFE = list_of(3)(6)(7)(1)(2)(0)(5)(4);
    const oligo::Permutate ABCDEFGH_ABCDEFGH(blockLength, ABCDEFGH, ABCDEFGH);
    const oligo::Permutate ABCDEFGH_ADBCEFGH(blockLength, ABCDEFGH, ADBCEFGH);
    const oligo::Permutate ABCDEFGH_DGHBCAFE(blockLength, ABCDEFGH, DGHBCAFE);
    const oligo::Permutate ADBCEFGH_DGHBCAFE(blockLength, ADBCEFGH, DGHBCAFE);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ABCDEFGH(kmer), kmer);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ADBCEFGH(kmer), 0xFE98DCBA76543210UL);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_DGHBCAFE(kmer), 0x983210DCBAFE5476UL);
    CPPUNIT_ASSERT_EQUAL(ADBCEFGH_DGHBCAFE(kmer), 0xDC3210BA98FE5476UL);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ABCDEFGH(abcdefgh), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ADBCEFGH(abcdefgh), adbcefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_DGHBCAFE(abcdefgh), dghbcafe);
    CPPUNIT_ASSERT_EQUAL(ADBCEFGH_DGHBCAFE(adbcefgh), dghbcafe);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ABCDEFGH.reorder(abcdefgh), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ADBCEFGH.reorder(adbcefgh), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_DGHBCAFE.reorder(dghbcafe), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ADBCEFGH_DGHBCAFE.reorder(dghbcafe), abcdefgh);
}

template <typename KmerT>
void testPermutate(KmerT original, const KmerT expected, const std::vector<isaac::oligo::Permutate> permutateList)
{
    CPPUNIT_ASSERT_EQUAL(original, permutateList.front()(original));
    CPPUNIT_ASSERT_EQUAL(original, permutateList.front().reorder(original));
    KmerT permuted = original;
    BOOST_FOREACH(const isaac::oligo::Permutate &permutate, permutateList)
    {
        permuted = permutate(permuted);
        CPPUNIT_ASSERT_EQUAL(original, permutate.reorder(permuted));
    }
    CPPUNIT_ASSERT_EQUAL(expected, permuted);
    CPPUNIT_ASSERT_EQUAL(original, permutateList.back().reorder(permuted));
}

const isaac::oligo::ShortKmerType ORIGINAL16(0x76543210U);
const isaac::oligo::ShortKmerType EXPECTED16(0x32107654U);
const isaac::oligo::KmerType ORIGINAL(0xFEDCBA9876543210UL);
const isaac::oligo::KmerType EXPECTED(0x76543210FEDCBA98UL);
const isaac::oligo::LongKmerType ORIGINAL64(isaac::oligo::LongKmerType(0x1111222233334444UL) << 64 | isaac::oligo::LongKmerType(0x5555666677778888UL));
const isaac::oligo::LongKmerType EXPECTED64(isaac::oligo::LongKmerType(0x5555666677778888UL) << 64 | isaac::oligo::LongKmerType(0x1111222233334444UL));
// This is needed for cases when isaac::oligo::Kmer is defined as __uint128_t or else CPPUNIT_ASSERT_EQUAL fails to compile
inline std::ostream & operator <<(std::ostream &os, const isaac::oligo::LongKmerType &kmer)
{
    return os << isaac::oligo::bases(kmer);
}


void TestPermutate::testTwoErrors()
{
    using namespace isaac;
    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::ShortKmerType>(2);
        CPPUNIT_ASSERT_EQUAL(6UL, permutateList.size());
        testPermutate(ORIGINAL16, EXPECTED16, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::KmerType>(2);
        CPPUNIT_ASSERT_EQUAL(6UL, permutateList.size());
        testPermutate(ORIGINAL, EXPECTED, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::LongKmerType>(2);
        CPPUNIT_ASSERT_EQUAL(6UL, permutateList.size());
        testPermutate(ORIGINAL64, EXPECTED64, permutateList);
    }
}

void TestPermutate::testFourErrors()
{
    using namespace isaac;
    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::ShortKmerType>(4);
        CPPUNIT_ASSERT_EQUAL(70UL, permutateList.size());
        testPermutate(ORIGINAL16, EXPECTED16, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::KmerType>(4);
        CPPUNIT_ASSERT_EQUAL(70UL, permutateList.size());
        testPermutate(ORIGINAL, EXPECTED, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::LongKmerType>(4);
        CPPUNIT_ASSERT_EQUAL(70UL, permutateList.size());
        testPermutate(ORIGINAL64, EXPECTED64, permutateList);
    }
}
