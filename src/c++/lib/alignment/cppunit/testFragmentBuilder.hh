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

#ifndef iSAAC_ALIGNMENT_TEST_FRAGMENT_BUILDER_HH
#define iSAAC_ALIGNMENT_TEST_FRAGMENT_BUILDER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <vector>

#include "alignment/FragmentBuilder.hh"
#include "flowcell/ReadMetadata.hh"
#include "flowcell/Layout.hh"

class TestFragmentBuilder : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestFragmentBuilder );
    CPPUNIT_TEST( testEmptyMatchList );
    CPPUNIT_TEST( testSingleSeed );
    CPPUNIT_TEST( testSeedOffset );
    CPPUNIT_TEST( testMultiSeed );
    CPPUNIT_TEST( testRepeats );
    CPPUNIT_TEST( testMismatches );
    CPPUNIT_TEST( testLeadingSoftClips );
    CPPUNIT_TEST( testTrailingSoftClips );
    CPPUNIT_TEST( testLeadingAndTrailingSoftClips );
    CPPUNIT_TEST_SUITE_END();
private:
    const isaac::flowcell::ReadMetadataList readMetadataList;
    const isaac::alignment::SeedMetadataList seedMetadataList;
    const isaac::flowcell::FlowcellLayoutList flowcells;
    const std::vector<isaac::reference::Contig> contigList;
    const std::vector<char> bcl0;
    const std::vector<char> bcl2;
    const std::vector<char> bcl3;
    const std::vector<char> bcl4l;
    const std::vector<char> bcl4t;
    const std::vector<char> bcl4lt;
    const unsigned tile0;
    const unsigned tile2;
    const unsigned clusterId0;
    const unsigned clusterId2;
    isaac::alignment::Cluster cluster0;
    isaac::alignment::Cluster cluster2;
    isaac::alignment::Cluster cluster3;
    isaac::alignment::Cluster cluster4l;
    isaac::alignment::Cluster cluster4t;
    isaac::alignment::Cluster cluster4lt;
    std::vector<isaac::alignment::Match> matchList;
    // auxiliary method for tests dedicated to a single seed on each read
    void auxSingleSeed(const unsigned s0, const unsigned s1);
public:
    TestFragmentBuilder();
    void setUp();
    void tearDown();
    void testEmptyMatchList();
    void testSingleSeed();
    void testSeedOffset();
    void testMultiSeed();
    void testRepeats();
    void testMismatches();
    void testLeadingSoftClips();
    void testTrailingSoftClips();
    void testLeadingAndTrailingSoftClips();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_FRAGMENT_BUILDER_HH
