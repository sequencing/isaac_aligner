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
