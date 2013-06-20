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
 ** \file testSampleSheetCsvGrammar.hh
 **
 ** tests boost::spirit grammar for SampleSheetCsv parsing
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_TEST_SAMPLE_SHEET_CSV_GRAMMAR_HH
#define iSAAC_OPTIONS_TEST_SAMPLE_SHEET_CSV_GRAMMAR_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

class TestSampleSheetCsvGrammar : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestSampleSheetCsvGrammar );
    CPPUNIT_TEST( testStandard );
    CPPUNIT_TEST( testDualBarcodeMiSeq );
    CPPUNIT_TEST( testAnotherDualBarcodeMiSeq);
    CPPUNIT_TEST( testNonMultiplexedMiSeq );
    CPPUNIT_TEST( testSingleBarcodeMiSeq );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testStandard();
    void testDualBarcodeMiSeq();
    void testAnotherDualBarcodeMiSeq();
    void testNonMultiplexedMiSeq();
    void testSingleBarcodeMiSeq();

};

#endif // #ifndef iSAAC_OPTIONS_TEST_SAMPLE_SHEET_CSV_GRAMMAR_HH

