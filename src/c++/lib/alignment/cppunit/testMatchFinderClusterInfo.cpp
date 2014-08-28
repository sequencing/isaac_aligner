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
#include "testMatchFinderClusterInfo.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestMatchFinderClusterInfo, registryName("MatchFinderClusterInfo"));

void TestMatchFinderClusterInfo::setUp()
{
}

void TestMatchFinderClusterInfo::tearDown()
{
}


void TestMatchFinderClusterInfo::testFields()
{
    using isaac::alignment::matchFinder::ClusterInfo;
    ClusterInfo none;
    ClusterInfo all;
    all.markReadComplete(0);
    all.markReadComplete(1);
    ClusterInfo other;
    other.setBarcodeIndex(0U);
    other.markReadComplete(1);
    ClusterInfo barcode;
    ClusterInfo r1;
    r1.markReadComplete(0);
    ClusterInfo r2;
    r2.markReadComplete(1);
    // barcode
    CPPUNIT_ASSERT(!barcode.isBarcodeSet());
    CPPUNIT_ASSERT_EQUAL(false, barcode.isReadComplete(0));
    CPPUNIT_ASSERT_EQUAL(false, barcode.isReadComplete(1));
    // r1
    CPPUNIT_ASSERT(!r1.isBarcodeSet());
    CPPUNIT_ASSERT_EQUAL(true, r1.isReadComplete(0));
    CPPUNIT_ASSERT_EQUAL(false, r1.isReadComplete(1));
    // r2
    CPPUNIT_ASSERT(!r2.isBarcodeSet());
    CPPUNIT_ASSERT_EQUAL(false, r2.isReadComplete(0));
    CPPUNIT_ASSERT_EQUAL(true, r2.isReadComplete(1));
    // none
    CPPUNIT_ASSERT(!none.isBarcodeSet());
    CPPUNIT_ASSERT_EQUAL(false, none.isReadComplete(0));
    CPPUNIT_ASSERT_EQUAL(false, none.isReadComplete(1));
    // all
    CPPUNIT_ASSERT(!all.isBarcodeSet());
    CPPUNIT_ASSERT_EQUAL(true, all.isReadComplete(0));
    CPPUNIT_ASSERT_EQUAL(true, all.isReadComplete(1));
    // OTHER
    CPPUNIT_ASSERT_EQUAL(0U, other.getBarcodeIndex());
    CPPUNIT_ASSERT(other.isBarcodeSet());
    CPPUNIT_ASSERT_EQUAL(false, other.isReadComplete(0));
    CPPUNIT_ASSERT_EQUAL(true, other.isReadComplete(1));

    other.setBarcodeIndex(1234U);
    CPPUNIT_ASSERT_EQUAL(1234U, other.getBarcodeIndex());
    CPPUNIT_ASSERT(other.isBarcodeSet());
    CPPUNIT_ASSERT_EQUAL(false, other.isReadComplete(0));
    CPPUNIT_ASSERT_EQUAL(true, other.isReadComplete(1));

    other.setBarcodeIndex(482U);
    CPPUNIT_ASSERT_EQUAL(482U, other.getBarcodeIndex());
    CPPUNIT_ASSERT(other.isBarcodeSet());
    CPPUNIT_ASSERT_EQUAL(false, other.isReadComplete(0));
    CPPUNIT_ASSERT_EQUAL(true, other.isReadComplete(1));

}

