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
 **/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>

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
    const oligo::Kmer kmer = 0xFEDCBA9876543210UL;
    const oligo::Kmer abcd = 0xFEDCBA9876543210UL;
    const oligo::Kmer adbc = 0xFEDC3210BA987654UL;
    const oligo::Kmer dbca = 0x3210BA987654FEDCUL;
    const unsigned blockLength = 8;
    const std::vector<unsigned> ABCD = list_of(0)(1)(2)(3);
    const std::vector<unsigned> ADBC = list_of(0)(3)(1)(2);
    const std::vector<unsigned> DBCA = list_of(3)(1)(2)(0);
    const oligo::Permutate ABCD_ABCD(blockLength, ABCD, ABCD);
    const oligo::Permutate ABCD_ADBC(blockLength, ABCD, ADBC);
    const oligo::Permutate ABCD_DBCA(blockLength, ABCD, DBCA);
    const oligo::Permutate ADBC_DBCA(blockLength, ADBC, DBCA);
    CPPUNIT_ASSERT_EQUAL(ABCD_ABCD(kmer), kmer);
    CPPUNIT_ASSERT_EQUAL(ABCD_ADBC(kmer), oligo::Kmer(0xFEDC3210BA987654UL));
    CPPUNIT_ASSERT_EQUAL(ABCD_DBCA(kmer), oligo::Kmer(0x3210BA987654FEDCUL));
    CPPUNIT_ASSERT_EQUAL(ADBC_DBCA(kmer), oligo::Kmer(0xBA9876543210FEDCUL));
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
    const oligo::Kmer kmer = 0xFEDCBA9876543210UL;
    const oligo::Kmer abcdefgh = 0xFEDCBA9876543210UL;
    const oligo::Kmer adbcefgh = 0xFE98DCBA76543210UL;
    const oligo::Kmer dghbcafe = 0x983210DCBAFE5476UL;
    const unsigned blockLength = 4;
    const std::vector<unsigned> ABCDEFGH = list_of(0)(1)(2)(3)(4)(5)(6)(7);
    const std::vector<unsigned> ADBCEFGH = list_of(0)(3)(1)(2)(4)(5)(6)(7);
    const std::vector<unsigned> DGHBCAFE = list_of(3)(6)(7)(1)(2)(0)(5)(4);
    const oligo::Permutate ABCDEFGH_ABCDEFGH(blockLength, ABCDEFGH, ABCDEFGH);
    const oligo::Permutate ABCDEFGH_ADBCEFGH(blockLength, ABCDEFGH, ADBCEFGH);
    const oligo::Permutate ABCDEFGH_DGHBCAFE(blockLength, ABCDEFGH, DGHBCAFE);
    const oligo::Permutate ADBCEFGH_DGHBCAFE(blockLength, ADBCEFGH, DGHBCAFE);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ABCDEFGH(kmer), kmer);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ADBCEFGH(kmer), oligo::Kmer(0xFE98DCBA76543210UL));
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_DGHBCAFE(kmer), oligo::Kmer(0x983210DCBAFE5476UL));
    CPPUNIT_ASSERT_EQUAL(ADBCEFGH_DGHBCAFE(kmer), oligo::Kmer(0xDC3210BA98FE5476UL));
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ABCDEFGH(abcdefgh), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ADBCEFGH(abcdefgh), adbcefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_DGHBCAFE(abcdefgh), dghbcafe);
    CPPUNIT_ASSERT_EQUAL(ADBCEFGH_DGHBCAFE(adbcefgh), dghbcafe);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ABCDEFGH.reorder(abcdefgh), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ADBCEFGH.reorder(adbcefgh), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_DGHBCAFE.reorder(dghbcafe), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ADBCEFGH_DGHBCAFE.reorder(dghbcafe), abcdefgh);
}

void TestPermutate::testTwoErrors()
{
    using namespace isaac;
    const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList(2);
    CPPUNIT_ASSERT_EQUAL(6UL, permutateList.size());
    const oligo::Kmer kmer(0xFEDCBA9876543210UL);
    CPPUNIT_ASSERT_EQUAL(kmer, permutateList.front()(kmer));
    CPPUNIT_ASSERT_EQUAL(kmer, permutateList.front().reorder(kmer));
    oligo::Kmer permuted = kmer;
    BOOST_FOREACH(const oligo::Permutate &permutate, permutateList)
    {
        permuted = permutate(permuted);
        CPPUNIT_ASSERT_EQUAL(kmer, permutate.reorder(permuted));
    }
    CPPUNIT_ASSERT_EQUAL(0x76543210FEDCBA98UL, permuted);
    CPPUNIT_ASSERT_EQUAL(kmer, permutateList.back().reorder(permuted));
}

void TestPermutate::testFourErrors()
{
    using namespace isaac;
    const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList(4);
    CPPUNIT_ASSERT_EQUAL(70UL, permutateList.size());
    const oligo::Kmer kmer(0xFEDCBA9876543210UL);
    CPPUNIT_ASSERT_EQUAL(kmer, permutateList.front()(kmer));
    CPPUNIT_ASSERT_EQUAL(kmer, permutateList.front().reorder(kmer));
    oligo::Kmer permuted = kmer;
    BOOST_FOREACH(const oligo::Permutate &permutate, permutateList)
    {
        permuted = permutate(permuted);
        CPPUNIT_ASSERT_EQUAL(kmer, permutate.reorder(permuted));
    }
    CPPUNIT_ASSERT_EQUAL(0x76543210FEDCBA98UL, permuted);
    CPPUNIT_ASSERT_EQUAL(kmer, permutateList.back().reorder(permuted));
}
