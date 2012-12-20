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
 **
 ** \file testSequencingAdapterListGrammar.cpp
 **
 ** tests boost::spirit grammar for use bases mask parsing
 **
 ** \author Roman Petrovski
 **/

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testSequencingAdapterListGrammar.hh"

#include "flowcell/SequencingAdapterListGrammar.hpp"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestSequencingAdapterListGrammar, registryName("SequencingAdapterListGrammar"));

void TestSequencingAdapterListGrammar::setUp()
{
}

void TestSequencingAdapterListGrammar::tearDown()
{
}

void testParsing(const std::string &test, const isaac::flowcell::SequencingAdapterMetadataList& expectedResult)
{
    std::string::const_iterator parseIt(test.begin());
    std::string::const_iterator parseEnd(test.end());

    isaac::flowcell::SequencingAdapterListGrammar<std::string::const_iterator> parser;
    isaac::flowcell::SequencingAdapterMetadataList result;
    if (!boost::spirit::qi::parse(parseIt, parseEnd, parser, result) || parseEnd != parseIt)
    {
        CPPUNIT_FAIL("Could not parse stream text:\n" + std::string(parseIt, parseEnd));
    }

    CPPUNIT_ASSERT_EQUAL(expectedResult.size(), result.size());
    for (unsigned i = 0; i < expectedResult.size(); ++i)
    {
        if (expectedResult.at(i) != result.at(i))
        {
            CPPUNIT_FAIL((boost::format("Parsed result %s of expression '%s' at position %d does not match the expected: %s") %
                result.at(i) % test % i % expectedResult.at(i)).str());
        }
    }
}

using isaac::flowcell::SequencingAdapterListGrammar;
typedef SequencingAdapterListGrammar<std::string::const_iterator> Parser;
void TestSequencingAdapterListGrammar::testNextera()
{
    std::string test("CTGTCTCTTATACACATCT*,*AGATGTGTATAAGAGACAG");
    const isaac::flowcell::SequencingAdapterMetadata testResult0("CTGTCTCTTATACACATCT", false, 0);
    const isaac::flowcell::SequencingAdapterMetadata testResult1("AGATGTGTATAAGAGACAG", true, 0);

    testParsing(test, boost::assign::list_of(testResult0)(testResult1));
    testParsing("Nextera", boost::assign::list_of(testResult0)(testResult1));
}

void TestSequencingAdapterListGrammar::testNexteraMp()
{
    std::string test("CTGTCTCTTATACACATCT,AGATGTGTATAAGAGACAG");
    const isaac::flowcell::SequencingAdapterMetadata testResult0("CTGTCTCTTATACACATCT", false);
    const isaac::flowcell::SequencingAdapterMetadata testResult1("AGATGTGTATAAGAGACAG", false);

    testParsing(test, boost::assign::list_of(testResult0)(testResult1));
    testParsing("NexteraMp", boost::assign::list_of(testResult0)(testResult1));
}

