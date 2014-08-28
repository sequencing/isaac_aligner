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

#ifndef iSAAC_ALIGNMENT_TEST_DUPLICATE_FILTERING_HH
#define iSAAC_ALIGNMENT_TEST_DUPLICATE_FILTERING_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <vector>

#include "build/FragmentIndex.hh"

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
    isaac::build::FStrandFragmentIndex fLeft1Frp_;
    isaac::build::FStrandFragmentIndex fLeft2Frp_;
    isaac::build::FStrandFragmentIndex fLeft3Frp_;
    isaac::build::RStrandOrShadowFragmentIndex rRight1Frp_;
    isaac::build::RStrandOrShadowFragmentIndex rRight2Frp_;
    isaac::build::RStrandOrShadowFragmentIndex rRight3Frp_;

    isaac::build::FStrandFragmentIndex fLeft1Ffp_;
    isaac::build::FStrandFragmentIndex fLeft2Ffp_;
    isaac::build::FStrandFragmentIndex fRight1Ffp_;
    isaac::build::FStrandFragmentIndex fRight2Ffp_;

    isaac::build::RStrandOrShadowFragmentIndex rLeft1Rrp_;
    isaac::build::RStrandOrShadowFragmentIndex rLeft2Rrp_;
    isaac::build::RStrandOrShadowFragmentIndex rRight1Rrp_;
    isaac::build::RStrandOrShadowFragmentIndex rRight2Rrp_;

    isaac::build::RStrandOrShadowFragmentIndex rLeft1Rfp_;
    isaac::build::RStrandOrShadowFragmentIndex rLeft2Rfp_;
    isaac::build::RStrandOrShadowFragmentIndex rLeft3Rfp_;
    isaac::build::FStrandFragmentIndex fRight1Rfp_;
    isaac::build::FStrandFragmentIndex fRight2Rfp_;
    isaac::build::FStrandFragmentIndex fRight3Rfp_;

    isaac::build::FStrandFragmentIndex fLeft1FrpMb1_;
    isaac::build::FStrandFragmentIndex fLeft2FrpMb1_;
    isaac::build::FStrandFragmentIndex fLeft3FrpMb0_;
    isaac::build::FStrandFragmentIndex fLeft4FrpMb0_;
    isaac::build::RStrandOrShadowFragmentIndex rRight1FrpMb1_;
    isaac::build::RStrandOrShadowFragmentIndex rRight2FrpMb1_;
    isaac::build::RStrandOrShadowFragmentIndex rRight3FrpMb0_;
    isaac::build::RStrandOrShadowFragmentIndex rRight4FrpMb0_;

    isaac::build::FStrandFragmentIndex f1Fsh_;
    isaac::build::FStrandFragmentIndex f2Fsh_;
    isaac::build::FStrandFragmentIndex f3Fsh_;
    isaac::build::FStrandFragmentIndex f4Fsh_;
    isaac::build::RStrandOrShadowFragmentIndex sh1Fsh_;
    isaac::build::RStrandOrShadowFragmentIndex sh2Fsh_;
    isaac::build::RStrandOrShadowFragmentIndex sh3Fsh_;
    isaac::build::RStrandOrShadowFragmentIndex sh4Fsh_;

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

