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
 **
 ** \file testSequencingAdapterListGrammar.hh
 **
 ** tests boost::spirit grammar for SequencingAdapterList parsing
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_TEST_SEQUENCING_ADAPTER_LIST_GRAMMAR_HH
#define iSAAC_OPTIONS_TEST_SEQUENCING_ADAPTER_LIST_GRAMMAR_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

class TestSequencingAdapterListGrammar : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestSequencingAdapterListGrammar );
    CPPUNIT_TEST( testNextera );
    CPPUNIT_TEST( testNexteraMp );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testNextera();
    void testNexteraMp();
};

#endif // #ifndef iSAAC_OPTIONS_TEST_SEQUENCING_ADAPTER_LIST_GRAMMAR_HH

