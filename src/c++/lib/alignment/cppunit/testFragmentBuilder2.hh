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

#ifndef iSAAC_ALIGNMENT_TEST_SEQUENCING_ADAPTER_HH
#define iSAAC_ALIGNMENT_TEST_SEQUENCING_ADAPTER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "alignment/FragmentBuilder.hh"
#include "alignment/SeedMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "alignment/matchSelector/SequencingAdapter.hh"

class TestFragmentBuilder2 : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestFragmentBuilder2 );
    CPPUNIT_TEST( testEverything );
    CPPUNIT_TEST_SUITE_END();
private:
    const isaac::flowcell::ReadMetadataList readMetadataList;
    const isaac::alignment::SeedMetadataList seedMetadataList;
    const isaac::flowcell::FlowcellLayoutList flowcells;
    isaac::alignment::Cigar cigarBuffer_;

public:
    TestFragmentBuilder2();
    void setUp();
    void tearDown();
    void testEverything();
    void testMismatchCount();
    void testMismatchCycles();
    void testMismatchCyclesWithSoftClip();
    void testGapped();
    void testGappedWithNs();

private:
    void align(
        const std::string &read,
        const std::string &reference,
        const isaac::alignment::matchSelector::SequencingAdapterList &adapters,
        isaac::alignment::FragmentMetadata &fragmentMetadata,
        const bool gapped = false);
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_SEQUENCING_ADAPTER_HH

