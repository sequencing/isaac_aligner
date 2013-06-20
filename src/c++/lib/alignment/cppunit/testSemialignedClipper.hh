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

#ifndef iSAAC_ALIGNMENT_TEST_SEMIALIGNED_CLIPPER_HH
#define iSAAC_ALIGNMENT_TEST_SEMIALIGNED_CLIPPER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "alignment/fragmentBuilder/UngappedAligner.hh"
#include "alignment/SeedMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "alignment/matchSelector/SemialignedEndsClipper.hh"
#include "alignment/matchSelector/SequencingAdapter.hh"

class TestSemialignedClipper : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestSemialignedClipper );
    CPPUNIT_TEST( testEverything );
    CPPUNIT_TEST_SUITE_END();
private:
    const isaac::flowcell::ReadMetadataList readMetadataList;
    const isaac::alignment::SeedMetadataList seedMetadataList;
    const isaac::flowcell::FlowcellLayoutList flowcells;
    isaac::alignment::fragmentBuilder::UngappedAligner ungappedAligner_;
    isaac::alignment::Cigar cigarBuffer_;
    isaac::alignment::matchSelector::SemialignedEndsClipper clipper_;

public:
    TestSemialignedClipper();
    void setUp();
    void tearDown();
    void testEverything();
    void testLeftClipForward();
    void testRightClipForward();
    void testRightClipStartBeforeRef();
    void testLeftClipStartBeforeRef();

private:
    void align(
        const std::string &read,
        const std::string &reference,
        const isaac::alignment::matchSelector::SequencingAdapterList &adapters,
        isaac::alignment::FragmentMetadata &fragmentMetadata);
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_SEMIALIGNED_CLIPPER_HH

