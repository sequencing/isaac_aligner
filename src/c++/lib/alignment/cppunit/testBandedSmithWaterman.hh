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
 **/

#ifndef iSAAC_ALIGNMENT_TEST_BANDED_SMITH_WATERMAN_HH
#define iSAAC_ALIGNMENT_TEST_BANDED_SMITH_WATERMAN_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <vector>

#include "alignment/BandedSmithWaterman.hh"

class TestBandedSmithWaterman : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestBandedSmithWaterman );
    //CPPUNIT_TEST( testCustom );
    CPPUNIT_TEST( testUngapped );
    CPPUNIT_TEST( testSingleInsertion );
    CPPUNIT_TEST( testSingleDeletion );
    CPPUNIT_TEST( testMultipleIndels );
    CPPUNIT_TEST( testOverflow );
    CPPUNIT_TEST_SUITE_END();
private:
    const isaac::alignment::BandedSmithWaterman bsw;
    const std::string genome;
public:
    TestBandedSmithWaterman();
    void setUp();
    void tearDown();
    void testCustom();
    void testUngapped();
    void testSingleInsertion();
    void testSingleDeletion();
    void testMultipleIndels();
    void testOverflow();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_BANDED_SMITH_WATERMAN_HH
