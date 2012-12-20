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
#include <boost/foreach.hpp>
#include <boost/assign.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testNeighborsFinder.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestNeighborsFinder, registryName("NeighborsFinder"));

void TestNeighborsFinder::setUp()
{
}

void TestNeighborsFinder::tearDown()
{
}

void TestNeighborsFinder::testFindNeighbors()
{
    using namespace isaac;
    using isaac::reference::NeighborsFinder;
    typedef NeighborsFinder::AnnotatedKmer AnnotatedKmer;
    typedef NeighborsFinder::KmerList KmerList;
    using boost::assign::list_of;
    const oligo::Kmer k0 = 0x1234567800000000UL;
    const oligo::Kmer k1 = 0x1234667800000000UL;
    std::vector<oligo::Kmer> kmerListAux = list_of
        (k0)
        (k0|0x55500000UL)
        (k0|0x55900000UL)
        (k0|0x55DAAA00UL)
        (k0|0x55A00000UL)
        (k0|0x55E00000UL)
        (k0|0x55F00000UL)
        (k0|0xFFFFFFFFUL)
        (k1)
        (k1|0x55500000UL)
        (k1|0x55900000UL)
        (k1|0x55D00000UL)
        (k1|0x55A00000UL)
        (k1|0x55EAAA00UL)
        (k1|0x55F00000UL)
        (k1|0xFFFFFFFFUL)
        ;
    KmerList kmerList;
    BOOST_FOREACH(oligo::Kmer kmer, kmerListAux)
    {
        kmerList.push_back(AnnotatedKmer(kmer, 0));
    }
    const KmerList::const_iterator i0 = kmerList.begin();
    const KmerList::const_iterator i1 = kmerList.begin() + 8;
    CPPUNIT_ASSERT_EQUAL(k0, i0->value);
    CPPUNIT_ASSERT_EQUAL(k1, i1->value);
    NeighborsFinder::findNeighbors(kmerList, 3);

    CPPUNIT_ASSERT_EQUAL(false, kmerList[0].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(true, kmerList[1].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(true, kmerList[2].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(false, kmerList[3].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(true, kmerList[4].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(true, kmerList[5].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(true, kmerList[6].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(false, kmerList[7].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(false, kmerList[8].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(true, kmerList[9].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(true, kmerList[10].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(true, kmerList[11].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(true, kmerList[12].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(false, kmerList[13].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(true, kmerList[14].hasNeighbors);
    CPPUNIT_ASSERT_EQUAL(false, kmerList[15].hasNeighbors);
}
