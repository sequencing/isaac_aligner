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


#include "RegistryName.hh"
#include "testMD5Sum.hh"

using namespace std;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestMD5Sum, registryName("MD5Sum"));

void TestMD5Sum::setUp()
{
}

void TestMD5Sum::tearDown()
{
}

// TODO, may want to move this into the test class or into MD5Sum
string calcMd5Sum(const std::string &input)
{
    // Process input
    isaac::common::MD5Sum md5sum;
    md5sum.update(input.c_str(), input.size());
    // Get the digest
    const isaac::common::MD5Sum::Digest digest = md5sum.getDigest();
    const string r = isaac::common::MD5Sum::toHexString(digest.data, 16);

    CPPUNIT_ASSERT_EQUAL( r.size(), std::string::size_type(32) );
    return r;
}

void TestMD5Sum::testMD5Digest()
{
    CPPUNIT_ASSERT_EQUAL( string("e59ff97941044f85df5297e1c302d260"), 
        calcMd5Sum("Hello World\n") );
    CPPUNIT_ASSERT_EQUAL( string("7fc56270e7a70fa81a5935b72eacbe29"), 
        calcMd5Sum("A") );
    CPPUNIT_ASSERT_EQUAL( string("36d5c9c8b25197c8c16e4e2d5ed9451e"), 
        calcMd5Sum("This is a longer test WITH ! some ^& weird characters \\\n") );
    // more tests here
}

