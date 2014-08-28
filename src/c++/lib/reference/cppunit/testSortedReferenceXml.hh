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

#ifndef iSAAC_REFERENCE_TEST_SORTED_REFERENCE_XML_HH
#define iSAAC_REFERENCE_TEST_SORTED_REFERENCE_XML_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "reference/SortedReferenceXml.hh"

class TestSortedReferenceXml : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestSortedReferenceXml );
    CPPUNIT_TEST( testAll );
    CPPUNIT_TEST( testWriter );
    CPPUNIT_TEST( testContigsOnly );
    CPPUNIT_TEST( testMasksOnly );
    CPPUNIT_TEST( testMerge );
    CPPUNIT_TEST_SUITE_END();
private:
    const std::string xmlString;
public:
    TestSortedReferenceXml();
    void setUp();
    void tearDown();
    void testHasNeighbors();
    void testAll();
    void testWriter();
    void testContigsOnly();
    void testMasksOnly();
    void testMerge();

    void checkContent(const isaac::reference::SortedReferenceMetadata &sortedReferenceMetadata);
    void checkContigs(const isaac::reference::SortedReferenceMetadata &sortedReferenceMetadata);
    void checkMasks(const isaac::reference::SortedReferenceMetadata &sortedReferenceMetadata);
};

#endif // #ifndef iSAAC_REFERENCE_TEST_SORTED_REFERENCE_XML_HH
