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
 **/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testSeedId.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestSeedId, registryName("SeedId"));

void TestSeedId::setUp()
{
}

void TestSeedId::tearDown()
{
}


void TestSeedId::testFields()
{
    using isaac::alignment::SeedId;
    const SeedId none(0, 0, 0, 0, 0);
    const SeedId all(SeedId::TILE_MASK, SeedId::BARCODE_MASK, SeedId::CLUSTER_MASK, SeedId::SEED_MASK, SeedId::REVERSE_MASK);
    SeedId other(4020, 1234UL, 1234567UL, 3, 1);
    const SeedId tile(SeedId::TILE_MASK, 0, 0, 0, 0);
    const SeedId barcode(0, SeedId::BARCODE_MASK, 0, 0, 0);
    const SeedId cluster(0, 0, SeedId::CLUSTER_MASK, 0, 0);
    const SeedId seed(0, 0, 0, SeedId::SEED_MASK, 0);
    const SeedId reverse(0, 0, 0, 0, SeedId::REVERSE_MASK);
    // tile
    CPPUNIT_ASSERT_EQUAL((unsigned long)SeedId::TILE_MASK, tile.getTile());
    CPPUNIT_ASSERT_EQUAL(0UL, tile.getCluster());
    CPPUNIT_ASSERT_EQUAL(0UL, tile.getReverse());
    CPPUNIT_ASSERT(!tile.isReverse());
    CPPUNIT_ASSERT_EQUAL(0UL, tile.getSeed());
    // barcode
    CPPUNIT_ASSERT_EQUAL((unsigned long)SeedId::BARCODE_MASK, barcode.getBarcode());
    // cluster
    CPPUNIT_ASSERT_EQUAL((unsigned long)SeedId::CLUSTER_MASK, cluster.getCluster());
    // reverse
    CPPUNIT_ASSERT_EQUAL((unsigned long)SeedId::REVERSE_MASK, reverse.getReverse());
    CPPUNIT_ASSERT_EQUAL(0UL, reverse.getTile());
    CPPUNIT_ASSERT_EQUAL(0UL, reverse.getCluster());
    CPPUNIT_ASSERT(reverse.isReverse());
    CPPUNIT_ASSERT_EQUAL(0UL, reverse.getSeed());
    // seed
    CPPUNIT_ASSERT_EQUAL((unsigned long)SeedId::SEED_MASK, seed.getSeed());
    CPPUNIT_ASSERT_EQUAL(0UL, seed.getTile());
    CPPUNIT_ASSERT_EQUAL(0UL, seed.getCluster());
    CPPUNIT_ASSERT_EQUAL(0UL, seed.getReverse());
    CPPUNIT_ASSERT(!seed.isReverse());
    // none
    CPPUNIT_ASSERT_EQUAL(0UL, none.getTile());
    CPPUNIT_ASSERT_EQUAL(0UL, none.getCluster());
    CPPUNIT_ASSERT_EQUAL(0UL, none.getReverse());
    CPPUNIT_ASSERT(!none.isReverse());
    CPPUNIT_ASSERT_EQUAL(0UL, none.getSeed());
    // all
    CPPUNIT_ASSERT_EQUAL((unsigned long)SeedId::TILE_MASK, all.getTile());
    CPPUNIT_ASSERT_EQUAL((unsigned long)SeedId::BARCODE_MASK, all.getBarcode());
    CPPUNIT_ASSERT_EQUAL((unsigned long)SeedId::CLUSTER_MASK, all.getCluster());
    CPPUNIT_ASSERT_EQUAL((unsigned long)SeedId::REVERSE_MASK, all.getReverse());
    CPPUNIT_ASSERT(all.isReverse());
    CPPUNIT_ASSERT_EQUAL((unsigned long)SeedId::SEED_MASK, all.getSeed());
    CPPUNIT_ASSERT_EQUAL(true, all.isNSeedId());
    // OTHER
    CPPUNIT_ASSERT_EQUAL(4020UL, other.getTile());
    CPPUNIT_ASSERT_EQUAL(1234UL, other.getBarcode());
    CPPUNIT_ASSERT_EQUAL(1234567UL, other.getCluster());
    CPPUNIT_ASSERT_EQUAL(1UL, other.getReverse());
    CPPUNIT_ASSERT(other.isReverse());
    CPPUNIT_ASSERT_EQUAL(3UL, other.getSeed());
    CPPUNIT_ASSERT_EQUAL(false, other.isNSeedId());
    // turning non-N-seed id into an N-seed id
    other.setNSeedId(false);
    CPPUNIT_ASSERT_EQUAL(4020UL, other.getTile());
    CPPUNIT_ASSERT_EQUAL(1234UL, other.getBarcode());
    CPPUNIT_ASSERT_EQUAL(1234567UL, other.getCluster());
    CPPUNIT_ASSERT_EQUAL(1UL, other.getReverse());
    CPPUNIT_ASSERT(other.isReverse());
    CPPUNIT_ASSERT_EQUAL((unsigned long)SeedId::SEED_MASK, other.getSeed());
    CPPUNIT_ASSERT_EQUAL(true, other.isNSeedId());
    // testing 'lowest seed' indicator
    other.setNSeedId(true);
    CPPUNIT_ASSERT_EQUAL(4020UL, other.getTile());
    CPPUNIT_ASSERT_EQUAL(1234UL, other.getBarcode());
    CPPUNIT_ASSERT_EQUAL(1234567UL, other.getCluster());
    CPPUNIT_ASSERT_EQUAL(0UL, other.getReverse());
    CPPUNIT_ASSERT(!other.isReverse());
    CPPUNIT_ASSERT_EQUAL((unsigned long)SeedId::SEED_MASK, other.getSeed());
    CPPUNIT_ASSERT_EQUAL(true, other.isNSeedId());
}

void TestSeedId::testOverflow()
{
    using isaac::alignment::SeedId;
    CPPUNIT_ASSERT_THROW(SeedId(SeedId::TILE_MASK + 1, 0, 0, 0, 0), isaac::common::PreConditionException);
    CPPUNIT_ASSERT_THROW(SeedId(0, SeedId::BARCODE_MASK + 1, 0, 0, 0), isaac::common::PreConditionException);
    CPPUNIT_ASSERT_THROW(SeedId(0, 0, SeedId::CLUSTER_MASK + 1, 0, 0), isaac::common::PreConditionException);
    CPPUNIT_ASSERT_THROW(SeedId(0, 0, 0, SeedId::SEED_MASK + 1, 0), isaac::common::PreConditionException);
    CPPUNIT_ASSERT_THROW(SeedId(0, 0, 0, 0, SeedId::REVERSE_MASK + 1), isaac::common::PreConditionException);
}

void TestSeedId::testSort()
{
    using isaac::alignment::SeedId;
    using std::vector;
    using std::sort;
    using boost::assign::list_of;
    std::vector<SeedId> v = list_of(SeedId(4,0,12,5,0))(SeedId(5,0,2,5,0))(SeedId(4,0,11,5,1))(SeedId(4,0,11,5,0));
    sort(v.begin(), v.end());
    CPPUNIT_ASSERT_EQUAL(4UL, v[0].getTile());
    CPPUNIT_ASSERT_EQUAL(4UL, v[1].getTile());
    CPPUNIT_ASSERT_EQUAL(4UL, v[2].getTile());
    CPPUNIT_ASSERT_EQUAL(5UL, v[3].getTile());
    CPPUNIT_ASSERT_EQUAL(11UL, v[0].getCluster());
    CPPUNIT_ASSERT_EQUAL(11UL, v[1].getCluster());
    CPPUNIT_ASSERT_EQUAL(12UL, v[2].getCluster());
    CPPUNIT_ASSERT(!v[0].isReverse());
    CPPUNIT_ASSERT(v[1].isReverse());
    CPPUNIT_ASSERT(!v[2].isReverse());
    CPPUNIT_ASSERT(!v[3].isReverse());
}
