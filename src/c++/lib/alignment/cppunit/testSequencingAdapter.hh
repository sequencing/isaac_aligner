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

#ifndef iSAAC_ALIGNMENT_TEST_SEQUENCING_ADAPTER_HH
#define iSAAC_ALIGNMENT_TEST_SEQUENCING_ADAPTER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "alignment/fragmentBuilder/UngappedAligner.hh"
#include "alignment/SeedMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "alignment/matchSelector/SequencingAdapter.hh"

class TestSequencingAdapter : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestSequencingAdapter );
    CPPUNIT_TEST( testEverything );
    CPPUNIT_TEST_SUITE_END();
private:
    const isaac::flowcell::ReadMetadataList readMetadataList;
    const isaac::alignment::SeedMetadataList seedMetadataList;
    const isaac::flowcell::FlowcellLayoutList flowcells;
    isaac::alignment::fragmentBuilder::UngappedAligner ungappedAligner_;
    const std::string irrelevantQualities;
    isaac::alignment::Cigar cigarBuffer_;

    isaac::alignment::matchSelector::SequencingAdapterList matePairAdapters;

    isaac::alignment::matchSelector::SequencingAdapterList standardAdapters;

public:
    TestSequencingAdapter();
    void setUp();
    void tearDown();
    void testEverything();
    void testMp51M49S();
    void testMp51S49M();
    void testMp94M6S();
    void testMp33S67M();
    void testMp40M60S();
    void testMp47S53M();
    void testMp30S70M();
    void testMp11S89M();
    void testMp16S84M();
    void testStd38M62S();
    void testStd76S24MReverse();
    void testStd36M114S();
    void testStdBeforeSequence();
    void testStdReverseAfterSequence();
    void testStdReverseSequenceTooGood();


private:
    void align(
        const std::string &read,
        const std::string &reference,
        const isaac::alignment::matchSelector::SequencingAdapterList &adapters,
        isaac::alignment::FragmentMetadata &fragmentMetadata);
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_SEQUENCING_ADAPTER_HH

