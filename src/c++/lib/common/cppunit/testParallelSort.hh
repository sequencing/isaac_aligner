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
 ** \file testParallelSort.hh
 **
 ** Unit tests for ParallelSort.hpp
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_COMMON_CPPUNIT_TEST_PARALLEL_SORT
#define iSAAC_COMMON_CPPUNIT_TEST_PARALLEL_SORT

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "common/ParallelSort.hpp"

class TestParallelSort : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestParallelSort );
    CPPUNIT_TEST( testSort );
    CPPUNIT_TEST_SUITE_END();
private:
    static const unsigned int BUFFER_SIZE = 1024;
    char buffer_[BUFFER_SIZE+1];
public:
    void setUp();
    void tearDown();
    void testSort();
};

#endif // #ifndef iSAAC_COMMON_CPPUNIT_TEST_PARALLEL_SORT
