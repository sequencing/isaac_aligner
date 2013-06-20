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

#ifndef iSAAC_COMMON_TEST_EXCEPTIONS_HH
#define iSAAC_COMMON_TEST_EXCEPTIONS_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "common/Exceptions.hh"

class TestExceptions : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestExceptions );
    CPPUNIT_TEST( testErrorNumber );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testErrorNumber();
};

#endif // #ifndef iSAAC_COMMON_TEST_EXCEPTIONS_HH
