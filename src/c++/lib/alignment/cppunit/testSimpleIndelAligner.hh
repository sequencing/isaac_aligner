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
 ** \file testSimpleIndelAligner.hh
 **
 ** Tests for medium-size gap aligner.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_TEST_SIMPLE_INDEL_ALIGNER_HH
#define iSAAC_ALIGNMENT_TEST_SIMPLE_INDEL_ALIGNER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "alignment/Cigar.hh"
#include "alignment/FragmentMetadata.hh"
#include "alignment/SeedMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "alignment/matchSelector/SequencingAdapter.hh"

class TestSimpleIndelAligner : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestSimpleIndelAligner );
    CPPUNIT_TEST( testEverything );
    CPPUNIT_TEST_SUITE_END();
private:
    isaac::alignment::Cigar cigarBuffer_;

public:
    TestSimpleIndelAligner();
    void setUp();
    void tearDown();
    void testEverything();

private:
    void align(
        const std::string &read,
        const std::string &reference,
        isaac::alignment::FragmentMetadataList &fragmentMetadataList);

    void align(
        const std::string &read,
        const std::string &reference,
        const isaac::alignment::SeedMetadataList &seedMetadataList,
        isaac::alignment::FragmentMetadataList &fragmentMetadataList);
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_SIMPLE_INDEL_ALIGNER_HH

