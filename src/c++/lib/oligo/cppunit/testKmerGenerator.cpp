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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>
#include <boost/utility/binary.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testKmerGenerator.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestKmerGenerator, registryName("KmerGenerator"));

void TestKmerGenerator::setUp()
{
}

void TestKmerGenerator::tearDown()
{
}


void TestKmerGenerator::testUnsigned()
{
    typedef isaac::oligo::KmerGenerator<unsigned> KmerGenerator;
    unsigned kmer;
    std::vector<char>::const_iterator position;
    {
        const std::string s("ANAACGTA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end(), 7);
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
    {
        const std::string s = std::string("ANAACGTAA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end(), 7);
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01B0U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 2L);
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
    {
        const std::string s = std::string("ANAACGTAANAAAAAACGT");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end(), 7);
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01B0U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 2L);
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 10L);
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x06U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 11L);
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01BU);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 12L);
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
}


void TestKmerGenerator::testConstMethods()
{
    typedef isaac::oligo::KmerGenerator<unsigned, std::string::const_iterator> KmerGenerator;
    unsigned kmer;
    std::vector<char>::const_iterator position;
    {
        const std::string s("AGAACGTA");
        KmerGenerator kmerGenerator(s.begin(), s.end(), 9);
        CPPUNIT_ASSERT_EQUAL(262143UL, isaac::oligo::getMaxKmer<unsigned long>(9));
        CPPUNIT_ASSERT(!isaac::oligo::generateKmer(9, kmer, s.begin(), s.end()));
    }
    {
        const std::string s = std::string("AACGTAA");
        KmerGenerator kmerGenerator(s.begin(), s.end(), 7);
        CPPUNIT_ASSERT_EQUAL(16383UL, isaac::oligo::getMaxKmer<unsigned long>(7));
        CPPUNIT_ASSERT(isaac::oligo::generateKmer(7, kmer, s.begin(), s.end()));
        CPPUNIT_ASSERT_EQUAL(kmer, unsigned(BOOST_BINARY(00 00 01 10 11 00 00)));
    }
}
