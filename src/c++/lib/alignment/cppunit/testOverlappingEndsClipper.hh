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

#ifndef iSAAC_ALIGNMENT_TEST_OVERLAPPING_ENDS_CLIPPER_HH
#define iSAAC_ALIGNMENT_TEST_OVERLAPPING_ENDS_CLIPPER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "alignment/matchSelector/OverlappingEndsClipper.hh"

class TestOverlappingEndsClipper : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestOverlappingEndsClipper );
    CPPUNIT_TEST( testEverything );
    CPPUNIT_TEST_SUITE_END();
private:

public:
    TestOverlappingEndsClipper();
    void setUp();
    void tearDown();
    void testEverything();

private:
    isaac::alignment::Cigar cigarBuffer_;
    isaac::alignment::Cluster cluster_;

    void init(
        const std::string &read1, const std::string &quality1, const bool read1Reverse,
        const std::string &read2, const std::string &quality2, const bool read2Reverse,
        const std::string &reference,
        isaac::alignment::BamTemplate &templ,
        isaac::reference::ContigList &contigList);
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_OVERLAPPING_ENDS_CLIPPER_HH

