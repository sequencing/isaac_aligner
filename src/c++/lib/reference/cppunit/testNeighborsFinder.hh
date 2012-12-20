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
 **/

#ifndef iSAAC_REFERENCE_TEST_NEIGHBORS_FINDER_HH
#define iSAAC_REFERENCE_TEST_NEIGHBORS_FINDER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "reference/NeighborsFinder.hh"

class TestNeighborsFinder : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestNeighborsFinder );
    CPPUNIT_TEST( testFindNeighbors );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testFindNeighbors();
};

#endif // #ifndef iSAAC_REFERENCE_TEST_NEIGHBORS_FINDER_HH
