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

#ifndef iSAAC_ALIGNMENT_TEST_SEED_ID_HH
#define iSAAC_ALIGNMENT_TEST_SEED_ID_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "alignment/SeedId.hh"

class TestSeedId : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestSeedId );
    CPPUNIT_TEST( testFields );
    CPPUNIT_TEST( testOverflow );
    CPPUNIT_TEST( testSort );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testFields();
    void testOverflow();
    void testSort();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_SEED_ID_HH

