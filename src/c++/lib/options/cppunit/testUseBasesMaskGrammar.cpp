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
 ** \file testUseBasesMaskGrammar.cpp
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
#include <boost/lexical_cast.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testUseBasesMaskGrammar.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestUseBasesMaskGrammar, registryName("UseBasesMaskGrammar"));

void TestUseBasesMaskGrammar::setUp()
{
}

void TestUseBasesMaskGrammar::tearDown()
{
}


using isaac::options::UseBasesMaskGrammar;
typedef UseBasesMaskGrammar<std::string::const_iterator> Parser;
void TestUseBasesMaskGrammar::testValid()
{
    std::vector<unsigned int> readLengths = boost::assign::list_of(10)(20)(30);

    {
        std::vector<std::string > result;
        std::string useBasesMask("y*n,y*n,y*n,Y*N");
        std::string::const_iterator parseIt(useBasesMask.begin());
        std::string::const_iterator parseEnd(useBasesMask.end());
        CPPUNIT_ASSERT(boost::spirit::qi::parse(parseIt, parseEnd, Parser(readLengths), result));
        CPPUNIT_ASSERT_MESSAGE("Could not parse the use-bases-mask: " + useBasesMask +
                                   " at pos:" + boost::lexical_cast<std::string>(parseIt - useBasesMask.begin()),
                               parseEnd == parseIt);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Expected '*' expansion to yield a mask of read length", readLengths[0], unsigned(result[0].size()));
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Expected '*' expansion to yield a mask of read length", readLengths[1], unsigned(result[1].size()));
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Expected '*' expansion to yield a mask of read length", readLengths[2], unsigned(result[2].size()));
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Expected code to assume read lenght 0 if it is not supplied", unsigned(2), unsigned(result[3].size()));
    }

    {
        std::vector<std::string > result;
        std::string useBasesMask("y10n,Y100n,y*n100");
        std::string::const_iterator parseIt(useBasesMask.begin());
        std::string::const_iterator parseEnd(useBasesMask.end());
        CPPUNIT_ASSERT(boost::spirit::qi::parse(parseIt, parseEnd, Parser(readLengths), result));
//        std::cerr << "use bases mask: " << useBasesMask << "\n";
//        std::cerr << "expanded use bases mask: " << boost::algorithm::join(result, ",") << "\n";
        CPPUNIT_ASSERT_MESSAGE("Could not parse the use-bases-mask: " + useBasesMask +
                                   " at pos:" + boost::lexical_cast<std::string>(parseIt - useBasesMask.begin()),
                               parseEnd == parseIt);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Expected non-* masks to ignore read length", unsigned(11), unsigned(result[0].size()));
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Expected non-* masks to ignore read length", unsigned(101), unsigned(result[1].size()));
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Expected '*' expansion to yield one mask base if the rest of the expression consumes entire read length",
                                     unsigned(101), unsigned(result[2].size()));
    }
}

void TestUseBasesMaskGrammar::testInvalid()
{
    std::vector<std::string > result;

    {
        std::vector<unsigned int> readLengths = boost::assign::list_of(10)(20)(30);
        std::string useBasesMask("y*k");
        std::string::const_iterator parseIt(useBasesMask.begin());
        std::string::const_iterator parseEnd(useBasesMask.end());
        boost::spirit::qi::parse(parseIt, parseEnd, Parser(readLengths), result);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Expected a parsing failure in use-bases-mask due to invalid character: " +
                                     useBasesMask + " at pos: 2", 2, int(parseIt - useBasesMask.begin()));
    }

    {
        std::vector<unsigned int> readLengths = boost::assign::list_of(10)(20)(30);
        std::string useBasesMask("y-1");
        std::string::const_iterator parseIt(useBasesMask.begin());
        std::string::const_iterator parseEnd(useBasesMask.end());
        boost::spirit::qi::parse(parseIt, parseEnd, Parser(readLengths), result);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Expected a parsing failure in use-bases-mask due to negative repeat count: " +
                                     useBasesMask + " at pos: 1", 1, int(parseIt - useBasesMask.begin()));
    }

    {
        std::vector<unsigned int> readLengths = boost::assign::list_of(10)(20)(30);
        std::string useBasesMask("y0");
        std::string::const_iterator parseIt(useBasesMask.begin());
        std::string::const_iterator parseEnd(useBasesMask.end());
        boost::spirit::qi::parse(parseIt, parseEnd, Parser(readLengths), result);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Expected a parsing failure in use-bases-mask due to 0 repeat count: " +
                                     useBasesMask + " at pos: 1", 1, int(parseIt - useBasesMask.begin()));
    }

    //    using isaac::alignment::SeedId;
//    CPPUNIT_ASSERT_THROW(SeedId(SeedId::TILE_MASK + 1, 0, 0, 0, 0), isaac::common::PreConditionException);
//    CPPUNIT_ASSERT_THROW(SeedId(0, SeedId::CLUSTER_MASK + 1, 0, 0, 0), isaac::common::PreConditionException);
//    CPPUNIT_ASSERT_THROW(SeedId(0, 0, SeedId::READ_MASK + 1, 0, 0), isaac::common::PreConditionException);
//    CPPUNIT_ASSERT_THROW(SeedId(0, 0, 0, SeedId::SEED_MASK + 1, 0), isaac::common::PreConditionException);
//    CPPUNIT_ASSERT_THROW(SeedId(0, 0, 0, 0, SeedId::REVERSE_MASK + 1), isaac::common::PreConditionException);
}

