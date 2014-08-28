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

#ifndef iSAAC_ALIGNMENT_TEST_GAP_REALIGNER_HH
#define iSAAC_ALIGNMENT_TEST_GAP_REALIGNER_HH

#include <cppunit/extensions/HelperMacros.h>

class TestGapRealigner : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestGapRealigner );
    CPPUNIT_TEST( testFull );
    CPPUNIT_TEST_SUITE_END();
private:

public:
    TestGapRealigner();
    void setUp();
    void tearDown();

    void testFull();
    void testMore();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_GAP_REALIGNER_HH

