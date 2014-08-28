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

#ifndef iSAAC_COMMON_TEST_MD5_HH
#define iSAAC_COMMON_TEST_MD5_HH

#include <cppunit/extensions/HelperMacros.h>
#include "common/MD5Sum.hh"

class TestMD5Sum : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestMD5Sum );
    CPPUNIT_TEST( testMD5Digest );
    CPPUNIT_TEST_SUITE_END();
public:
    void setUp();
    void tearDown();
    void testMD5Digest();
};

#endif // #ifndef iSAAC_COMMON_TEST_MD5_HH
