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
 **
 ** \file testBarcodeResolver.hh
 **
 ** tests internals of barcode resolver
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_TEST_BARCODE_RESOLVER_HH
#define iSAAC_OPTIONS_TEST_BARCODE_RESOLVER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

class TestBarcodeResolver : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestBarcodeResolver );
    CPPUNIT_TEST( testGet1MismatchKmer );
    CPPUNIT_TEST( testOneComponent );
    CPPUNIT_TEST( testTwoComponents );
    CPPUNIT_TEST( testMismatchCollision );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testGet1MismatchKmer();
    void testOneComponent();
    void testTwoComponents();
    void testMismatchCollision();
};

#endif // #ifndef iSAAC_OPTIONS_TEST_BARCODE_RESOLVER_HH

