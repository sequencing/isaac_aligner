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

    void testAllTogether();

    void testFull();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_GAP_REALIGNER_HH

