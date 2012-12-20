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
 ** \file testFastIo.cpp
 **
 ** Unit tests for FastIo.hh
 **
 ** \author Come Raczy
 **/

#include <iostream>
#include <sstream>
#include <string>

using namespace std;

#include "RegistryName.hh"
#include "testFastIo.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestFastIo, registryName("FastIo"));

void TestFastIo::setUp()
{
    memset(buffer_, 0, BUFFER_SIZE+1);
}

void TestFastIo::tearDown()
{
}


void TestFastIo::testSprintFloatZeros()
{
    char *buffer = buffer_;
    isaac::common::sprintFloat<0, 0, 10>(buffer, 0.0);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0"));
    isaac::common::sprintFloat<1, 0, 10>(buffer, 0.0);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.0"));
    isaac::common::sprintFloat<2, 0, 10>(buffer, 0.0);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.00"));
    isaac::common::sprintFloat<5, 0, 10>(buffer, 0.0);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.00000"));
}

void TestFastIo::testSprintFloatZerosPadding()
{
    char *buffer = buffer_;
    isaac::common::sprintFloat<0, 2, 10>(buffer, -0.0);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string(" 0"));
    *buffer = 0;
    isaac::common::sprintFloat<5, 6, 10>(buffer, -0.0);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.00000"));
    *buffer = 0;
    isaac::common::sprintFloat<5, 10, 10>(buffer, -0.0);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("   0.00000"));
    *buffer = 0;
    isaac::common::sprintFloat<5, 10, 10>(buffer, -0.0);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("   0.00000"));
    *buffer = 0;
}

void TestFastIo::testSprintFloatSmallPositive()
{
    char *buffer = buffer_;
    isaac::common::sprintFloat<0, 0, 10>(buffer, 0.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0"));
    *buffer = 0;
    isaac::common::sprintFloat<1, 0, 10>(buffer, 0.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.0"));
    *buffer = 0;
    isaac::common::sprintFloat<2, 0, 10>(buffer, 0.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.04"));
    *buffer = 0;
    isaac::common::sprintFloat<3, 0, 10>(buffer, 0.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.040"));
    *buffer = 0;
    isaac::common::sprintFloat<3, 0, 10>(buffer, 0.049450001);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.049"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, 0.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.0400"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, 0.049);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.0490"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, 0.0494);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.0494"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, 0.04944);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.0494"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, 0.049450001);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.0495"));
    *buffer = 0;
}

void TestFastIo::testSprintFloatSmallNegative()
{
    char *buffer = buffer_;
    isaac::common::sprintFloat<0, 0, 10>(buffer, -0.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0"));
    *buffer = 0;
    isaac::common::sprintFloat<1, 0, 10>(buffer, -0.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("0.0"));
    *buffer = 0;
    isaac::common::sprintFloat<2, 0, 10>(buffer, -0.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-0.04"));
    *buffer = 0;
    isaac::common::sprintFloat<3, 0, 10>(buffer, -0.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-0.040"));
    *buffer = 0;
    isaac::common::sprintFloat<3, 0, 10>(buffer, -0.049450001);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-0.049"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, -0.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-0.0400"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, -0.049);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-0.0490"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, -0.0494);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-0.0494"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, -0.04944);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-0.0494"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, -0.049450001);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-0.0495"));
    *buffer = 0;
}

void TestFastIo::testSprintFloatMediumPositive()
{
    char *buffer = buffer_;
    isaac::common::sprintFloat<0, 0, 10>(buffer, 6.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("6"));
    *buffer = 0;
    isaac::common::sprintFloat<1, 0, 10>(buffer, 6.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("6.0"));
    *buffer = 0;
    isaac::common::sprintFloat<2, 0, 10>(buffer, 6.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("6.04"));
    *buffer = 0;
    isaac::common::sprintFloat<3, 0, 10>(buffer, 6.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("6.040"));
    *buffer = 0;
    isaac::common::sprintFloat<3, 0, 10>(buffer, 6.049450001);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("6.049"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, 6.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("6.0400"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, 6.049);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("6.0490"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, 6.0494);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("6.0494"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, 6.04944);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("6.0494"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, 4.049450001);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("4.0495"));
    *buffer = 0;
}

void TestFastIo::testSprintFloatMediumNegative()
{
    char *buffer = buffer_;
    isaac::common::sprintFloat<0, 0, 10>(buffer, -6.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-6"));
    *buffer = 0;
    isaac::common::sprintFloat<1, 0, 10>(buffer, -6.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-6.0"));
    *buffer = 0;
    isaac::common::sprintFloat<2, 0, 10>(buffer, -6.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-6.04"));
    *buffer = 0;
    isaac::common::sprintFloat<3, 0, 10>(buffer, -6.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-6.040"));
    *buffer = 0;
    isaac::common::sprintFloat<3, 0, 10>(buffer, -6.049450001);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-6.049"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, -6.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-6.0400"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, -6.049);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-6.0490"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, -6.0494);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-6.0494"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, -6.04944);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-6.0494"));
    *buffer = 0;
    isaac::common::sprintFloat<4, 0, 10>(buffer, -1.049450001);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-1.0495"));
    *buffer = 0;
}

void TestFastIo::testSprintFloatLargePositive()
{
    char *buffer = buffer_;
    isaac::common::sprintFloat<0, 0, 10>(buffer, 1789012.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("1789012"));
    *buffer = 0;
    isaac::common::sprintFloat<1, 0, 10>(buffer, 1789012.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("1789012.0"));
    *buffer = 0;
    isaac::common::sprintFloat<2, 0, 10>(buffer, 1789012.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer,6), string("1789012.04",6));
    *buffer = 0;
    isaac::common::sprintFloat<3, 0, 10>(buffer, 1789012.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer,5), string("789012.040",5)); // truncated
    *buffer = 0;
    isaac::common::sprintFloat<3, 0, 10>(buffer, 1789012.049450001);
    CPPUNIT_ASSERT_EQUAL(string(buffer,5), string("789012.049",5));
    if (4 < sizeof(long))
    {
        *buffer = 0;
        isaac::common::sprintFloat<4, 0, 10>(buffer, 6789012.04);
        CPPUNIT_ASSERT_EQUAL(string(buffer,4), string("89012.0400",4));
        *buffer = 0;
        isaac::common::sprintFloat<4, 0, 10>(buffer, 6789012.049);
        CPPUNIT_ASSERT_EQUAL(string(buffer,4), string("89012.0490",4));
        *buffer = 0;
        isaac::common::sprintFloat<4, 0, 10>(buffer, 6789012.0494);
        CPPUNIT_ASSERT_EQUAL(string(buffer,4), string("89012.0494",4));
        *buffer = 0;
        isaac::common::sprintFloat<4, 0, 10>(buffer, 6789012.04944);
        CPPUNIT_ASSERT_EQUAL(string(buffer,4), string("89012.0494",4));
        *buffer = 0;
        isaac::common::sprintFloat<4, 0, 10>(buffer, 4789012.049450001);
        CPPUNIT_ASSERT_EQUAL(string(buffer,4), string("89012.0495",4));
        *buffer = 0;
    }
}

void TestFastIo::testSprintFloatLargeNegative()
{
    char *buffer = buffer_;
    isaac::common::sprintFloat<0, 0, 10>(buffer, -1789012.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-1789012"));
    *buffer = 0;
    isaac::common::sprintFloat<1, 0, 10>(buffer, -1789012.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer), string("-1789012.0"));
    *buffer = 0;
    isaac::common::sprintFloat<2, 0, 10>(buffer, -1789012.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer,6), string("1789012.04",6));
    *buffer = 0;
    isaac::common::sprintFloat<3, 0, 10>(buffer, -1789012.04);
    CPPUNIT_ASSERT_EQUAL(string(buffer,5), string("789012.040",5));
    *buffer = 0;
    isaac::common::sprintFloat<3, 0, 10>(buffer, -1789012.049450001);
    CPPUNIT_ASSERT_EQUAL(string(buffer,5), string("789012.049",5));
    if (4 < sizeof(long))
    {
        *buffer = 0;
        isaac::common::sprintFloat<4, 0, 10>(buffer, -6789012.04);
        CPPUNIT_ASSERT_EQUAL(string(buffer,4), string("89012.0400",4));
        *buffer = 0;
        isaac::common::sprintFloat<4, 0, 10>(buffer, -6789012.049);
        CPPUNIT_ASSERT_EQUAL(string(buffer,4), string("89012.0490",4));
        *buffer = 0;
        isaac::common::sprintFloat<4, 0, 10>(buffer, -6789012.0494);
        CPPUNIT_ASSERT_EQUAL(string(buffer,4), string("89012.0494",4));
        *buffer = 0;
        isaac::common::sprintFloat<4, 0, 10>(buffer, -6789012.04944);
        CPPUNIT_ASSERT_EQUAL(string(buffer,4), string("89012.0494",4));
        *buffer = 0;
        isaac::common::sprintFloat<4, 0, 10>(buffer, -1789012.049450001);
        CPPUNIT_ASSERT_EQUAL(string(buffer,4), string("89012.0495",4));
        *buffer = 0;
    }
}

void TestFastIo::testPutUnsignedInteger()
{
    std::ostringstream os;
    isaac::common::putUnsignedInteger<unsigned char>(os, 0);
    CPPUNIT_ASSERT_EQUAL(string("0"), os.str());
    os.str("");
    isaac::common::putUnsignedInteger<unsigned char>(os, 1);
    CPPUNIT_ASSERT_EQUAL(string("1"), os.str());
    os.str("");
    isaac::common::putUnsignedInteger<unsigned char>(os, 23);
    CPPUNIT_ASSERT_EQUAL(string("23"), os.str());
    os.str("");
    isaac::common::putUnsignedInteger<unsigned char>(os, 127);
    CPPUNIT_ASSERT_EQUAL(string("127"), os.str());
    os.str("");
    isaac::common::putUnsignedInteger<unsigned char>(os, 128);
    CPPUNIT_ASSERT_EQUAL(string("128"), os.str());
    os.str("");
    isaac::common::putUnsignedInteger<unsigned char>(os, 255);
    CPPUNIT_ASSERT_EQUAL(string("255"), os.str());
    os.str("");
    isaac::common::putUnsignedInteger<unsigned char>(os, static_cast<unsigned char>(256));
    CPPUNIT_ASSERT_EQUAL(string("0"), os.str());
    os.str("");
    isaac::common::putUnsignedInteger<unsigned char>(os, static_cast<unsigned char>(-1));
    CPPUNIT_ASSERT_EQUAL(string("255"), os.str());
    os.str("");
    // testing on unsigned short
    isaac::common::putUnsignedInteger<unsigned short>(os, 60000);
    CPPUNIT_ASSERT_EQUAL(string("60000"), os.str());
    os.str("");
    // testing on unsigned int
    isaac::common::putUnsignedInteger<unsigned int>(os, 123456789);
    CPPUNIT_ASSERT_EQUAL(string("123456789"), os.str());
    os.str("");
    // testing on unsigned default does not compile
    //isaac::common::putUnsignedInteger(os, 123456789);
    //CPPUNIT_ASSERT_EQUAL(string("123456789"), os.str());
    //os.str("");
}

void TestFastIo::testAppendUnsignedInteger()
{
    using isaac::common::appendUnsignedInteger;
    std::string expected;
    std::string s = expected;
    appendUnsignedInteger(s,0);
    expected.append("0");
    CPPUNIT_ASSERT_EQUAL(expected, s);
    appendUnsignedInteger(s,1);
    expected.append("1");
    CPPUNIT_ASSERT_EQUAL(expected, s);
    appendUnsignedInteger(s,10);
    expected.append("10");
    CPPUNIT_ASSERT_EQUAL(expected, s);
    appendUnsignedInteger(s,11);
    expected.append("11");
    CPPUNIT_ASSERT_EQUAL(expected, s);
    appendUnsignedInteger(s,0);
    expected.append("0");
    CPPUNIT_ASSERT_EQUAL(expected, s);
    appendUnsignedInteger(s,1234);
    expected.append("1234");
    CPPUNIT_ASSERT_EQUAL(expected, s);
}

void TestFastIo::testPutInteger()
{
    std::ostringstream os;
    isaac::common::putInteger<char>(os, 0);
    CPPUNIT_ASSERT_EQUAL(string("0"), os.str());
    os.str("");    
    isaac::common::putInteger<char>(os, -0);
    CPPUNIT_ASSERT_EQUAL(string("0"), os.str());
    os.str("");
    isaac::common::putInteger<char>(os, -1);
    CPPUNIT_ASSERT_EQUAL(string("-1"), os.str());
    os.str("");
    isaac::common::putInteger<char>(os, -128);
    CPPUNIT_ASSERT_EQUAL(string("-128"), os.str());
    os.str("");
    isaac::common::putInteger<char>(os, 127);
    CPPUNIT_ASSERT_EQUAL(string("127"), os.str());
    os.str("");
}

void TestFastIo::testGetUnsignedInteger()
{
    // tests for unsigned chars
    {
        std::istringstream is;
        typedef unsigned char Unsigned;
        Unsigned value=98;
        is.clear();
        is.str("0");
        isaac::common::getUnsignedInteger<Unsigned>(is, value);
        CPPUNIT_ASSERT_EQUAL(0, (int)value);
        is.clear();
        is.str("00000");
        isaac::common::getUnsignedInteger<Unsigned>(is, value);
        CPPUNIT_ASSERT_EQUAL(0, (int)value);
        is.clear();
        is.str("010");
        isaac::common::getUnsignedInteger<Unsigned>(is, value);
        CPPUNIT_ASSERT_EQUAL(10, (int)value);
        is.clear();
        is.str("1");
        isaac::common::getUnsignedInteger<Unsigned>(is, value);
        CPPUNIT_ASSERT_EQUAL(1, (int)value);
        is.clear();
        is.str("255");
        isaac::common::getUnsignedInteger<Unsigned>(is, value);
        CPPUNIT_ASSERT_EQUAL(255, (int)value);
        is.clear();
        is.str("256");
        isaac::common::getUnsignedInteger<Unsigned>(is, value);
        CPPUNIT_ASSERT_EQUAL(0, (int)value);
        is.clear();
        // this function does not read '-'
        value = 8;
        is.str("-1");
        isaac::common::getUnsignedInteger<Unsigned>(is, value);
        CPPUNIT_ASSERT_EQUAL(0, (int)value);
    }
    // tests for unsigned int
    {
        std::istringstream is;
        typedef unsigned int Unsigned;
        Unsigned value=98765;
        is.clear();
        is.str("0");
        isaac::common::getUnsignedInteger<Unsigned>(is, value);
        CPPUNIT_ASSERT_EQUAL(0, (int)value);
        is.clear();
        is.str("00000");
        isaac::common::getUnsignedInteger<Unsigned>(is, value);
        CPPUNIT_ASSERT_EQUAL(0, (int)value);
        is.clear();
        is.str("000012345000");
        isaac::common::getUnsignedInteger<Unsigned>(is, value);
        CPPUNIT_ASSERT_EQUAL(12345000, (int)value);
    }
}

void TestFastIo::testGetInteger()
{
    // tests for unsigned chars
    {
        std::istringstream is;
        typedef char Signed;
        Signed value=98;
        is.clear();
        is.str("-0001");
        isaac::common::getInteger<Signed>(is, value);
        CPPUNIT_ASSERT_EQUAL(-1, (int)value);
        is.clear();
        is.str("-128");
        isaac::common::getInteger<Signed>(is, value);
        CPPUNIT_ASSERT_EQUAL(-128, (int)value);
        is.clear();
        is.str("-129");
        isaac::common::getInteger<Signed>(is, value);
        CPPUNIT_ASSERT_EQUAL(127, (int)value);
        is.clear();
        is.str("128");
        isaac::common::getInteger<Signed>(is, value);
        CPPUNIT_ASSERT_EQUAL(-128, (int)value);
        is.clear();
        is.str("45abc");
        isaac::common::getInteger<Signed>(is, value);
        CPPUNIT_ASSERT_EQUAL(45, (int)value);
        char c = is.get();
        CPPUNIT_ASSERT_EQUAL('a', c);
        is.clear();
        is.str("45abc");
        isaac::common::getInteger<Signed>(is, value, true);
        CPPUNIT_ASSERT_EQUAL(45, (int)value);
        c = is.get();
        CPPUNIT_ASSERT_EQUAL('b', c);
    }
    // tests for unsigned int
    {
        std::istringstream is;
        typedef unsigned int Signed;
        Signed value=98765;
        is.clear();
        is.str("-0");
        isaac::common::getInteger<Signed>(is, value);
        CPPUNIT_ASSERT_EQUAL(0, (int)value);
        value = 765;
        is.clear();
        is.str("-00000");
        isaac::common::getInteger<Signed>(is, value);
        CPPUNIT_ASSERT_EQUAL(0, (int)value);
        is.clear();
        is.str("-123456");
        isaac::common::getInteger<Signed>(is, value);
        CPPUNIT_ASSERT_EQUAL(-123456, (int)value);
        // this method does not process nor choke on garbage
        is.clear();
        is.str("--123");
        value = 543;
        isaac::common::getInteger<Signed>(is, value);
        CPPUNIT_ASSERT_EQUAL(0, (int)value);
    }
}

void TestFastIo::testBoolIo()
{
    std::istringstream is("10YN ");
    bool b = false;
    isaac::common::getBool<'1', '0'>(is, b);
    CPPUNIT_ASSERT(b);
    isaac::common::getBool<'1', '0'>(is, b);
    CPPUNIT_ASSERT(!b);
    isaac::common::getBool<'Y', 'N'>(is, b);
    CPPUNIT_ASSERT(b);
    isaac::common::getBool<'Y', 'N'>(is, b);
    CPPUNIT_ASSERT(!b);
    CPPUNIT_ASSERT(is);
    isaac::common::getBool<'Y', 'N'>(is, b);
    CPPUNIT_ASSERT(is.fail());
    is.clear();
    is.str("");
    CPPUNIT_ASSERT(is);
    isaac::common::getBool<'Y', 'N'>(is, b);
    CPPUNIT_ASSERT(is.eof());
    // output
    std::ostringstream os;
    isaac::common::putBool<'1','0'>(os, true);
    isaac::common::putBool<'1','0'>(os, false);
    isaac::common::putBool<'Y','N'>(os, true);
    isaac::common::putBool<'Y','N'>(os, false);
    CPPUNIT_ASSERT_EQUAL(std::string("10YN"), os.str());
}
