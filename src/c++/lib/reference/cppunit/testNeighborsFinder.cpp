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

// This is needed for cases when isaac::oligo::Kmer is defined as __uint128_t or else CPPUNIT_ASSERT_EQUAL fails to compile
inline std::ostream & operator <<(std::ostream &os, const __uint128_t &kmer)
{
    return isaac::oligo::traceKmerValue(os, kmer);
}

template <typename KmerT>
void testFindNeighbors(const KmerT &MASK0, const KmerT &MASK1)
{
    using namespace isaac;
    using isaac::reference::NeighborsFinder;
    typedef typename NeighborsFinder<KmerT>::AnnotatedKmer AnnotatedKmer;
    typedef typename NeighborsFinder<KmerT>::KmerList KmerList;
    using boost::assign::list_of;
    std::vector<KmerT> kmerListAux = list_of
        (MASK0)
        (MASK0|0x55500000UL)
        (MASK0|0x55900000UL)
        (MASK0|0x55DAAA00UL)
        (MASK0|0x55A00000UL)
        (MASK0|0x55E00000UL)
        (MASK0|0x55F00000UL)
        (MASK0|0xFFFFFFFFUL)
        (MASK1)
        (MASK1|0x55500000UL)
        (MASK1|0x55900000UL)
        (MASK1|0x55D00000UL)
        (MASK1|0x55A00000UL)
        (MASK1|0x55EAAA00UL)
        (MASK1|0x55F00000UL)
        (MASK1|0xFFFFFFFFUL)
        ;
    KmerList kmerList;
    BOOST_FOREACH(KmerT kmer, kmerListAux)
    {
        kmerList.push_back(AnnotatedKmer(kmer, 0));
    }
    const typename KmerList::const_iterator i0 = kmerList.begin();
    const typename KmerList::const_iterator i1 = kmerList.begin() + 8;
    CPPUNIT_ASSERT_EQUAL(MASK0, i0->value);
    CPPUNIT_ASSERT_EQUAL(MASK1, i1->value);
    NeighborsFinder<KmerT>::findNeighbors(kmerList, 3);

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

void TestNeighborsFinder::testFindNeighbors()
{
    {
        const isaac::oligo::KmerType MASK0 = isaac::oligo::KmerType(0x1234567800000000UL);
        const isaac::oligo::KmerType MASK1 = isaac::oligo::KmerType(0x1234667800000000UL);
        ::testFindNeighbors(MASK0, MASK1);
    }

    {
        const isaac::oligo::LongKmerType MASK0 = isaac::oligo::LongKmerType(0x1234567800000000UL) << 64;
        const isaac::oligo::LongKmerType MASK1 = isaac::oligo::LongKmerType(0x1234667800000000UL) << 64;
        ::testFindNeighbors(MASK0, MASK1);
    }
}


