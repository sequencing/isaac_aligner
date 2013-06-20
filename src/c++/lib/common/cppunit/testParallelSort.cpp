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
 **
 ** \file testParallelSort.cpp
 **
 ** Unit tests for ParallelSort.hpp
 **
 ** \author Come Raczy
 **/

#include <iostream>
#include <sstream>
#include <string>

using namespace std;

#include "RegistryName.hh"
#include "testParallelSort.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestParallelSort, registryName("ParallelSort"));

void TestParallelSort::setUp()
{
    memset(buffer_, 0, BUFFER_SIZE+1);
}

void TestParallelSort::tearDown()
{
}


void TestParallelSort::testSort()
{
    std::vector<int> v;
    for (unsigned int i = 0; 101 > i; ++i)
    {
        v.push_back(rand());
    }
    std::vector<int> vv(v);
    std::sort(vv.begin(), vv.end(), std::less<int>());
    isaac::common::parallelSort(v, std::less<int>());
    CPPUNIT_ASSERT_EQUAL(v.size(), vv.size());
    for (size_t i = 0; v.size() > i; ++i)
    {
        CPPUNIT_ASSERT_EQUAL(v[i], vv[i]);
    }
}
