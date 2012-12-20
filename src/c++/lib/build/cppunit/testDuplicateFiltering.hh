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

#ifndef iSAAC_ALIGNMENT_TEST_DUPLICATE_FILTERING_HH
#define iSAAC_ALIGNMENT_TEST_DUPLICATE_FILTERING_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <vector>

#include "io/FragmentIndex.hh"

class TestDuplicateFiltering : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestDuplicateFiltering );
    CPPUNIT_TEST( testFrp );
    CPPUNIT_TEST( testFfp );
    CPPUNIT_TEST( testRrp );
    CPPUNIT_TEST( testRfp );
    CPPUNIT_TEST( testFrpReverseMatesInDifferentBins );
    CPPUNIT_TEST( testFsh );
    CPPUNIT_TEST( testAllTogether );
    CPPUNIT_TEST_SUITE_END();
private:
    isaac::io::FStrandFragmentIndex fLeft1Frp_;
    isaac::io::FStrandFragmentIndex fLeft2Frp_;
    isaac::io::FStrandFragmentIndex fLeft3Frp_;
    isaac::io::RStrandOrShadowFragmentIndex rRight1Frp_;
    isaac::io::RStrandOrShadowFragmentIndex rRight2Frp_;
    isaac::io::RStrandOrShadowFragmentIndex rRight3Frp_;

    isaac::io::FStrandFragmentIndex fLeft1Ffp_;
    isaac::io::FStrandFragmentIndex fLeft2Ffp_;
    isaac::io::FStrandFragmentIndex fRight1Ffp_;
    isaac::io::FStrandFragmentIndex fRight2Ffp_;

    isaac::io::RStrandOrShadowFragmentIndex rLeft1Rrp_;
    isaac::io::RStrandOrShadowFragmentIndex rLeft2Rrp_;
    isaac::io::RStrandOrShadowFragmentIndex rRight1Rrp_;
    isaac::io::RStrandOrShadowFragmentIndex rRight2Rrp_;

    isaac::io::RStrandOrShadowFragmentIndex rLeft1Rfp_;
    isaac::io::RStrandOrShadowFragmentIndex rLeft2Rfp_;
    isaac::io::RStrandOrShadowFragmentIndex rLeft3Rfp_;
    isaac::io::FStrandFragmentIndex fRight1Rfp_;
    isaac::io::FStrandFragmentIndex fRight2Rfp_;
    isaac::io::FStrandFragmentIndex fRight3Rfp_;

    isaac::io::FStrandFragmentIndex fLeft1FrpMb1_;
    isaac::io::FStrandFragmentIndex fLeft2FrpMb1_;
    isaac::io::FStrandFragmentIndex fLeft3FrpMb0_;
    isaac::io::FStrandFragmentIndex fLeft4FrpMb0_;
    isaac::io::RStrandOrShadowFragmentIndex rRight1FrpMb1_;
    isaac::io::RStrandOrShadowFragmentIndex rRight2FrpMb1_;
    isaac::io::RStrandOrShadowFragmentIndex rRight3FrpMb0_;
    isaac::io::RStrandOrShadowFragmentIndex rRight4FrpMb0_;

    isaac::io::FStrandFragmentIndex f1Fsh_;
    isaac::io::FStrandFragmentIndex f2Fsh_;
    isaac::io::FStrandFragmentIndex f3Fsh_;
    isaac::io::FStrandFragmentIndex f4Fsh_;
    isaac::io::RStrandOrShadowFragmentIndex sh1Fsh_;
    isaac::io::RStrandOrShadowFragmentIndex sh2Fsh_;
    isaac::io::RStrandOrShadowFragmentIndex sh3Fsh_;
    isaac::io::RStrandOrShadowFragmentIndex sh4Fsh_;

public:
    TestDuplicateFiltering();
    void setUp();
    void tearDown();
    void testFrp();
    void testFfp();
    void testRrp();
    void testRfp();
    void testFrpReverseMatesInDifferentBins();
    void testFsh();
    void testAllTogether();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_DUPLICATE_FILTERING_HH

