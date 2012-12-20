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

using namespace std;

#include "RegistryName.hh"
#include "testBarcodeId.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestBarcodeId, registryName("BarcodeId"));

void TestBarcodeId::setUp()
{
}

void TestBarcodeId::tearDown()
{
}


void TestBarcodeId::testFields()
{
    using isaac::demultiplexing::BarcodeId;
    const BarcodeId none(0, 0, 0, 0);
    const BarcodeId all(BarcodeId::TILE_MASK, BarcodeId::BARCODE_MASK, BarcodeId::CLUSTER_MASK, BarcodeId::MISMATCHES_MASK);
//    BarcodeId other(4020, 482UL, 1234567UL);
    BarcodeId other(123UL, 482UL, 0x2298a, 2);
    const BarcodeId tile(BarcodeId::TILE_MASK, 0, 0, 0);
    const BarcodeId barcode(0, BarcodeId::BARCODE_MASK, 0, 0);
    const BarcodeId cluster(0, 0, BarcodeId::CLUSTER_MASK, 0);
    const BarcodeId mismatches(0, 0, 0, BarcodeId::MISMATCHES_MASK);
    // tile
    CPPUNIT_ASSERT_EQUAL((unsigned long)BarcodeId::TILE_MASK, tile.getTile());
    CPPUNIT_ASSERT_EQUAL(0UL, tile.getBarcode());
    CPPUNIT_ASSERT_EQUAL(0UL, tile.getCluster());
    CPPUNIT_ASSERT_EQUAL(0UL, tile.getMismatches());
    // barcode
    CPPUNIT_ASSERT_EQUAL(0UL, barcode.getTile());
    CPPUNIT_ASSERT_EQUAL((unsigned long)BarcodeId::BARCODE_MASK, barcode.getBarcode());
    CPPUNIT_ASSERT_EQUAL(0UL, barcode.getCluster());
    CPPUNIT_ASSERT_EQUAL(0UL, barcode.getMismatches());
    // cluster
    CPPUNIT_ASSERT_EQUAL(0UL, cluster.getTile());
    CPPUNIT_ASSERT_EQUAL(0UL, cluster.getBarcode());
    CPPUNIT_ASSERT_EQUAL((unsigned long)BarcodeId::CLUSTER_MASK, cluster.getCluster());
    CPPUNIT_ASSERT_EQUAL(0UL, cluster.getMismatches());
    // mismatches
    CPPUNIT_ASSERT_EQUAL(0UL, mismatches.getTile());
    CPPUNIT_ASSERT_EQUAL(0UL, mismatches.getBarcode());
    CPPUNIT_ASSERT_EQUAL(0UL, mismatches.getCluster());
    CPPUNIT_ASSERT_EQUAL((unsigned long)BarcodeId::MISMATCHES_MASK, mismatches.getMismatches());
    // none
    CPPUNIT_ASSERT_EQUAL(0UL, none.getTile());
    CPPUNIT_ASSERT_EQUAL(0UL, none.getBarcode());
    CPPUNIT_ASSERT_EQUAL(0UL, none.getCluster());
    CPPUNIT_ASSERT_EQUAL(0UL, none.getCluster());
    // all
    CPPUNIT_ASSERT_EQUAL((unsigned long)BarcodeId::TILE_MASK, all.getTile());
    CPPUNIT_ASSERT_EQUAL((unsigned long)BarcodeId::BARCODE_MASK, all.getBarcode());
    CPPUNIT_ASSERT_EQUAL((unsigned long)BarcodeId::CLUSTER_MASK, all.getCluster());
    CPPUNIT_ASSERT_EQUAL((unsigned long)BarcodeId::MISMATCHES_MASK, all.getMismatches());
    // OTHER
    CPPUNIT_ASSERT_EQUAL(123UL, other.getTile());
    CPPUNIT_ASSERT_EQUAL(482UL, other.getBarcode());
    CPPUNIT_ASSERT_EQUAL(0x2298aUL, other.getCluster());
    CPPUNIT_ASSERT_EQUAL(2UL, other.getMismatches());
}

void TestBarcodeId::testOverflow()
{
    using isaac::demultiplexing::BarcodeId;
    CPPUNIT_ASSERT_THROW(BarcodeId(BarcodeId::TILE_MASK + 1, 0, 0, 0), isaac::common::PreConditionException);
    CPPUNIT_ASSERT_THROW(BarcodeId(0, BarcodeId::BARCODE_MASK + 1, 0, 0), isaac::common::PreConditionException);
    CPPUNIT_ASSERT_THROW(BarcodeId(0, 0, BarcodeId::CLUSTER_MASK + 1, 0), isaac::common::PreConditionException);
    CPPUNIT_ASSERT_THROW(BarcodeId(0, 0, 0, BarcodeId::MISMATCHES_MASK + 1), isaac::common::PreConditionException);
}
