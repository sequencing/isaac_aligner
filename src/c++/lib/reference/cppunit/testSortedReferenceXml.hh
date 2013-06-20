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
