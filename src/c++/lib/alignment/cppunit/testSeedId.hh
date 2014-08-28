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

