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

#ifndef iSAAC_ALIGNMENT_TEST_CLUSTER_INFO_HH
#define iSAAC_ALIGNMENT_TEST_CLUSTER_INFO_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "alignment/matchFinder/TileClusterInfo.hh"

class TestMatchFinderClusterInfo : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestMatchFinderClusterInfo );
    CPPUNIT_TEST( testFields );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testFields();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_CLUSTER_INFO_HH

