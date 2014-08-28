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
 ** \file testUseBasesMaskGrammar.hh
 **
 ** tests boost::spirit grammar for use bases mask parsing
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_TEST_USE_BASES_MASK_GRAMMAR_HH
#define iSAAC_OPTIONS_TEST_USE_BASES_MASK_GRAMMAR_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "options/UseBasesMaskGrammar.hh"

class TestUseBasesMaskGrammar : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestUseBasesMaskGrammar );
    CPPUNIT_TEST( testValid );
    CPPUNIT_TEST( testInvalid );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testValid();
    void testInvalid();
};

#endif // #ifndef iSAAC_OPTIONS_TEST_USE_BASES_MASK_GRAMMAR_HH

