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

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>
#include <boost/assign/std/vector.hpp> 

#include "build/DuplicatePairEndFilter.hh"
#include "build/DuplicateFragmentIndexFiltering.hh"

using namespace std;
using namespace isaac::io;
using namespace isaac::build;
using namespace isaac::reference;


#include "RegistryName.hh"
#include "testDuplicateFiltering.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestDuplicateFiltering, registryName("TestDuplicateFiltering"));


class FakePackedFragmentBuffer : public PackedFragmentBuffer
{
public:
    void fillWithUniqueClusterIdPattern()
    {
        for(std::vector<char>::iterator it = begin(); end() != it; it += sizeof(isaac::io::FragmentHeader))
        {
            isaac::io::FragmentHeader &header = *(isaac::io::FragmentHeader*)&*it;
            header.clusterId_ = std::distance(begin(), it);
        }
    }
};

template<typename IndexT> class TestDuplicateFilter
{

};

isaac::build::BarcodeBamMapping::BarcodeSampleIndexMap emptyMap;
template<>
class TestDuplicateFilter<isaac::build::RStrandOrShadowFragmentIndex> : public RSDuplicateFilter<false>
{
public:
    TestDuplicateFilter() : RSDuplicateFilter<false>(emptyMap){}
};

template<>
class TestDuplicateFilter<isaac::build::FStrandFragmentIndex> : public FDuplicateFilter<false>
{
public:
    TestDuplicateFilter() : FDuplicateFilter<false>(emptyMap){}
};

template <typename IndexT>
void testNoDifferences(std::vector<IndexT> bin, std::vector<IndexT>  expectedUnique)
{
    typedef IndexT IndexType;
    DuplicatePairEndFilter filter(false);

    isaac::alignment::BinMetadataList binMetadataList(1);
    isaac::alignment::BinMetadataCRefList binMetadataCRefList(1, boost::cref(binMetadataList.front()));
    binMetadataList[0] = isaac::alignment::BinMetadata(0, 0, isaac::reference::ReferencePosition(0,0), 1000, "", 0);
    // Resize to fit the highest offset fragment we have in the bin
    binMetadataList.at(0).incrementDataSize(isaac::reference::ReferencePosition(0,0), 10000 * sizeof(isaac::io::FragmentHeader));
    FakePackedFragmentBuffer fakeEmptyFragmentBuffer;
    fakeEmptyFragmentBuffer.resize(binMetadataList.at(0));
    fakeEmptyFragmentBuffer.fillWithUniqueClusterIdPattern();

    isaac::flowcell::BarcodeMetadataList barcodeMetadataList(1);
    BuildStats fakeBuildStats(binMetadataCRefList, barcodeMetadataList);
    std::vector<PackedFragmentBuffer::Index> filteredIndex;
    filter.filterInput(
        TestDuplicateFilter<IndexT>(), fakeEmptyFragmentBuffer, bin.begin(), bin.end(), fakeBuildStats, 0, std::back_inserter(filteredIndex));

    std::vector<unsigned long> uniqueFragments;
    std::transform(filteredIndex.begin(), filteredIndex.end(), std::back_inserter(uniqueFragments),
                   boost::bind(&PackedFragmentBuffer::Index::dataOffset_, _1));
    std::sort(uniqueFragments.begin(), uniqueFragments.end());

    std::vector<unsigned long> expectedUniqueFragments;
    std::transform(expectedUnique.begin(), expectedUnique.end(), std::back_inserter(expectedUniqueFragments),
                   boost::bind(&IndexType::dataOffset_, _1));
    std::sort(expectedUniqueFragments.begin(), expectedUniqueFragments.end());

    std::vector<unsigned long> diff;
    std::set_symmetric_difference(uniqueFragments.begin(), uniqueFragments.end(),
                        expectedUniqueFragments.begin(), expectedUniqueFragments.end(),
                        std::back_inserter(diff));
    if (!diff.empty())
    {
        std::cerr << "unexpected difference " << diff[0] << std::endl;
    }
    CPPUNIT_ASSERT_EQUAL(size_t(0), diff.size());

}
/**
 * \brief Set up the fragment pairs. Naming convention:
 * <strand><relative location><pair number><pair orientation>
 * fragments with pair number 1 are expected to be the best choices
 *
 * dataOffset_ field is used for unique identification during results comparison
 */
TestDuplicateFiltering::TestDuplicateFiltering():
        //common pairs
        fLeft1Frp_(ReferencePosition(0, 0),
                   FragmentIndexMate(false, true, 1, FragmentIndexAnchor(300)), 3),
        fLeft2Frp_(ReferencePosition(0, 0),
                   FragmentIndexMate(false, true, 1, FragmentIndexAnchor(300)), 2),
        fLeft3Frp_(ReferencePosition(0, 0),
                   FragmentIndexMate(false, true, 1, FragmentIndexAnchor(300)), 3),

        rRight1Frp_(ReferencePosition(0, 200), FragmentIndexAnchor(300),
                    FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 3),
        rRight2Frp_(ReferencePosition(0, 200), FragmentIndexAnchor(300),
                    FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 2),
        rRight3Frp_(ReferencePosition(0, 200), FragmentIndexAnchor(300),
                    FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 3),

        fLeft1Ffp_(ReferencePosition(0, 0),
                   FragmentIndexMate(false, false, 0, FragmentIndexAnchor(200)), 3),
        fLeft2Ffp_(ReferencePosition(0, 0),
                   FragmentIndexMate(false, false, 0, FragmentIndexAnchor(200)), 2),

        fRight1Ffp_(ReferencePosition(0, 200),
                    FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 3),
        fRight2Ffp_(ReferencePosition(0, 200),
                    FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 2),

        rLeft1Rrp_(ReferencePosition(0, 0), FragmentIndexAnchor(100),
                   FragmentIndexMate(false, true, 1, FragmentIndexAnchor(300)), 3),
        rLeft2Rrp_(ReferencePosition(0, 0), FragmentIndexAnchor(100),
                   FragmentIndexMate(false, true, 1, FragmentIndexAnchor(300)), 2),

        rRight1Rrp_(ReferencePosition(0, 200), FragmentIndexAnchor(300),
                    FragmentIndexMate(false, true, 1, FragmentIndexAnchor(100)), 3),
        rRight2Rrp_(ReferencePosition(0, 200), FragmentIndexAnchor(300),
                    FragmentIndexMate(false, true, 1, FragmentIndexAnchor(100)), 2),

        rLeft1Rfp_(ReferencePosition(0, 0), FragmentIndexAnchor(100),
                   FragmentIndexMate(false, false, 0, FragmentIndexAnchor(200)), 3),
        rLeft2Rfp_(ReferencePosition(0, 0), FragmentIndexAnchor(100),
                   FragmentIndexMate(false, false, 0, FragmentIndexAnchor(200)), 2),
        rLeft3Rfp_(ReferencePosition(0, 0), FragmentIndexAnchor(100),
                   FragmentIndexMate(false, false, 0, FragmentIndexAnchor(200)), 3),

        fRight1Rfp_(ReferencePosition(0, 200),
                    FragmentIndexMate(false, true, 1, FragmentIndexAnchor(100)), 3),
        fRight2Rfp_(ReferencePosition(0, 200),
                    FragmentIndexMate(false, true, 1, FragmentIndexAnchor(100)), 2),
        fRight3Rfp_(ReferencePosition(0, 200),
                    FragmentIndexMate(false, true, 1, FragmentIndexAnchor(100)), 3),

        // pairs with reverse-stranded mates stored in different bins
        fLeft1FrpMb1_(ReferencePosition(0, 0),
                      FragmentIndexMate(false, true, 1, FragmentIndexAnchor(300)), 3),
        fLeft2FrpMb1_(ReferencePosition(0, 0),
                      FragmentIndexMate(false, true, 1, FragmentIndexAnchor(300)), 2),
        fLeft3FrpMb0_(ReferencePosition(0, 0),
                      FragmentIndexMate(false, true, 0, FragmentIndexAnchor(300)), 3),
        fLeft4FrpMb0_(ReferencePosition(0, 0),
                      FragmentIndexMate(false, true, 0, FragmentIndexAnchor(300)), 3),

        rRight1FrpMb1_(ReferencePosition(0, 200), FragmentIndexAnchor(300),
                       FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 3),
        rRight2FrpMb1_(ReferencePosition(0, 200), FragmentIndexAnchor(300),
                       FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 2),
        rRight3FrpMb0_(ReferencePosition(0, 200), FragmentIndexAnchor(300),
                       FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 3),
        rRight4FrpMb0_(ReferencePosition(0, 200), FragmentIndexAnchor(300),
                       FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 3),

//        // FSingleton-shadow cases
        f1Fsh_(ReferencePosition(0, 0),
               FragmentIndexMate(true, false, 0, FragmentIndexAnchor(0x0123012301230123)), 3),
        f2Fsh_(ReferencePosition(0, 0),
               FragmentIndexMate(true, false, 0, FragmentIndexAnchor(0x0123012301230123)), 2),
        f3Fsh_(ReferencePosition(0, 0),
               FragmentIndexMate(true, false, 0, FragmentIndexAnchor(0x1230123012301230)), 3),
        f4Fsh_(ReferencePosition(0, 0),
               FragmentIndexMate(true, false, 0, FragmentIndexAnchor(0x1230123012301230)), 3),

        sh1Fsh_(ReferencePosition(0, 0), FragmentIndexAnchor(0x0123012301230123),
                FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 3),
        sh2Fsh_(ReferencePosition(0, 0), FragmentIndexAnchor(0x0123012301230123),
                FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 2),
        sh3Fsh_(ReferencePosition(0, 0), FragmentIndexAnchor(0x1230123012301230),
                FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 3),
        sh4Fsh_(ReferencePosition(0, 0), FragmentIndexAnchor(0x1230123012301230),
                FragmentIndexMate(false, false, 0, FragmentIndexAnchor(0)), 3)

{
    fLeft1Frp_.dataOffset_ = 0 * sizeof(isaac::io::FragmentHeader);
    fLeft2Frp_.dataOffset_ = 1 * sizeof(isaac::io::FragmentHeader);
    fLeft3Frp_.dataOffset_ = 2 * sizeof(isaac::io::FragmentHeader);

    rRight1Frp_.dataOffset_ = 1000 * sizeof(isaac::io::FragmentHeader);
    rRight2Frp_.dataOffset_ = 1001 * sizeof(isaac::io::FragmentHeader);
    rRight3Frp_.dataOffset_ = 1002 * sizeof(isaac::io::FragmentHeader);

    fLeft1Ffp_.dataOffset_ = 10 * sizeof(isaac::io::FragmentHeader);
    fLeft2Ffp_.dataOffset_ = 11 * sizeof(isaac::io::FragmentHeader);

    fRight1Ffp_.dataOffset_ = 1010 * sizeof(isaac::io::FragmentHeader);
    fRight2Ffp_.dataOffset_ = 1011 * sizeof(isaac::io::FragmentHeader);

    rLeft1Rrp_.dataOffset_ = 20 * sizeof(isaac::io::FragmentHeader);
    rLeft2Rrp_.dataOffset_ = 21 * sizeof(isaac::io::FragmentHeader);

    rRight1Rrp_.dataOffset_ = 1020 * sizeof(isaac::io::FragmentHeader);
    rRight2Rrp_.dataOffset_ = 1021 * sizeof(isaac::io::FragmentHeader);

    rLeft1Rfp_.dataOffset_ = 30 * sizeof(isaac::io::FragmentHeader);
    rLeft2Rfp_.dataOffset_ = 31 * sizeof(isaac::io::FragmentHeader);
    rLeft3Rfp_.dataOffset_ = 32 * sizeof(isaac::io::FragmentHeader);

    fRight1Rfp_.dataOffset_ = 1030 * sizeof(isaac::io::FragmentHeader);
    fRight2Rfp_.dataOffset_ = 1031 * sizeof(isaac::io::FragmentHeader);
    fRight3Rfp_.dataOffset_ = 1032 * sizeof(isaac::io::FragmentHeader);

    // pairs with reverse-stranded mates stored in different bins
    fLeft1FrpMb1_.dataOffset_ = 40 * sizeof(isaac::io::FragmentHeader);
    fLeft2FrpMb1_.dataOffset_ = 41 * sizeof(isaac::io::FragmentHeader);
    fLeft3FrpMb0_.dataOffset_ = 42 * sizeof(isaac::io::FragmentHeader);
    fLeft4FrpMb0_.dataOffset_ = 43 * sizeof(isaac::io::FragmentHeader);

    rRight1FrpMb1_.dataOffset_ = 1040 * sizeof(isaac::io::FragmentHeader);
    rRight2FrpMb1_.dataOffset_ = 1041 * sizeof(isaac::io::FragmentHeader);
    rRight3FrpMb0_.dataOffset_ = 1042 * sizeof(isaac::io::FragmentHeader);
    rRight4FrpMb0_.dataOffset_ = 1043 * sizeof(isaac::io::FragmentHeader);

//        // FSingleton-shadow cases
    f1Fsh_.dataOffset_ = 50 * sizeof(isaac::io::FragmentHeader);
    f2Fsh_.dataOffset_ = 51 * sizeof(isaac::io::FragmentHeader);
    f3Fsh_.dataOffset_ = 52 * sizeof(isaac::io::FragmentHeader);
    f4Fsh_.dataOffset_ = 53 * sizeof(isaac::io::FragmentHeader);

    sh1Fsh_.dataOffset_ = 1050 * sizeof(isaac::io::FragmentHeader);
    sh2Fsh_.dataOffset_ = 1051 * sizeof(isaac::io::FragmentHeader);
    sh3Fsh_.dataOffset_ = 1052 * sizeof(isaac::io::FragmentHeader);
    sh4Fsh_.dataOffset_ = 1053 * sizeof(isaac::io::FragmentHeader);

}

void TestDuplicateFiltering::setUp()
{
}

void TestDuplicateFiltering::tearDown()
{
}

void TestDuplicateFiltering::testFrp()
{
    std::vector<isaac::build::FStrandFragmentIndex> fInput = boost::assign::list_of
        (fLeft1Frp_)(fLeft2Frp_);
    std::vector<isaac::build::FStrandFragmentIndex> fExpectedResults = boost::assign::list_of
        (fLeft1Frp_);
    testNoDifferences(fInput, fExpectedResults);

    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsInput = boost::assign::list_of
        (rRight1Frp_)(rRight2Frp_);
    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsExpectedResults = boost::assign::list_of
        (rRight1Frp_);
    testNoDifferences(rsInput, rsExpectedResults);
}

void TestDuplicateFiltering::testFfp()
{
    std::vector<isaac::build::FStrandFragmentIndex> fInput = boost::assign::list_of
        (fLeft1Ffp_)(fLeft2Ffp_)
        (fRight1Ffp_)(fRight2Ffp_);
    std::vector<isaac::build::FStrandFragmentIndex> fExpectedResults = boost::assign::list_of
        (fLeft1Ffp_)
        (fRight1Ffp_);
    testNoDifferences(fInput, fExpectedResults);
}

void TestDuplicateFiltering::testRrp()
{
    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsInput = boost::assign::list_of
        (rLeft1Rrp_)(rLeft2Rrp_)
        (rRight1Rrp_)(rRight2Rrp_);

    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsExpectedResults = boost::assign::list_of
        (rLeft1Rrp_)
        (rRight1Rrp_);

    testNoDifferences(rsInput, rsExpectedResults);
}

void TestDuplicateFiltering::testRfp()
{
    std::vector<isaac::build::FStrandFragmentIndex> fInput = boost::assign::list_of
        (fRight1Rfp_)(fRight2Rfp_);
    std::vector<isaac::build::FStrandFragmentIndex> fExpectedResults = boost::assign::list_of
        (fRight1Rfp_);
    testNoDifferences(fInput, fExpectedResults);

    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsInput = boost::assign::list_of
        (rLeft1Rfp_)(rLeft2Rfp_);
    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsExpectedResults = boost::assign::list_of
        (rLeft1Rfp_);
    testNoDifferences(rsInput, rsExpectedResults);
}


void TestDuplicateFiltering::testFrpReverseMatesInDifferentBins()
{
    std::vector<isaac::build::FStrandFragmentIndex> fInput = boost::assign::list_of
        (fLeft1FrpMb1_)(fLeft2FrpMb1_)(fLeft3FrpMb0_)(fLeft4FrpMb0_);
    // left mates are in the same bin, 2 must surwive out of 4
    std::vector<isaac::build::FStrandFragmentIndex> fExpectedResults = boost::assign::list_of
        (fLeft1FrpMb1_)(fLeft3FrpMb0_);
    testNoDifferences(fInput, fExpectedResults);

    // the reverse mates are in different bins (by the test case definition)

    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsInput1 = boost::assign::list_of
        (rRight1FrpMb1_)(rRight2FrpMb1_);
    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsExpectedResults1 = boost::assign::list_of
        (rRight1FrpMb1_);
    testNoDifferences(rsInput1, rsExpectedResults1);

    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsInput0 = boost::assign::list_of
        (rRight3FrpMb0_)(rRight4FrpMb0_);
    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsExpectedResults0 = boost::assign::list_of
        (rRight3FrpMb0_);

    testNoDifferences(rsInput0, rsExpectedResults0);
}

void TestDuplicateFiltering::testFsh()
{
    std::vector<isaac::build::FStrandFragmentIndex> fInput = boost::assign::list_of
        (f1Fsh_)(f2Fsh_)(f3Fsh_)(f4Fsh_);
    // left mates are in the same bin, 2 must surwive out of 4
    std::vector<isaac::build::FStrandFragmentIndex> fExpectedResults = boost::assign::list_of
        (f1Fsh_)(f3Fsh_);
    testNoDifferences(fInput, fExpectedResults);

    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsInput = boost::assign::list_of
        (sh1Fsh_)(sh2Fsh_)(sh3Fsh_)(sh4Fsh_);
    // the reverse mates are in different bins (by the test case definition)
    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsExpectedResults = boost::assign::list_of
        (sh1Fsh_)(sh3Fsh_);
    testNoDifferences(rsInput, rsExpectedResults);

}

void TestDuplicateFiltering::testAllTogether()
{
    std::vector<isaac::build::FStrandFragmentIndex> fInput = boost::assign::list_of
        (fLeft1Frp_)(fLeft2Frp_)(fLeft3Frp_)
        (fLeft1Ffp_)(fLeft2Ffp_)
        (fRight1Ffp_)(fRight2Ffp_)
        (fRight1Rfp_)(fRight2Rfp_)(fRight3Rfp_)
        (fLeft1FrpMb1_)(fLeft2FrpMb1_)(fLeft3FrpMb0_)(fLeft4FrpMb0_)
        (f1Fsh_)(f2Fsh_)(f3Fsh_)(f4Fsh_);
    std::vector<isaac::build::FStrandFragmentIndex> fExpectedResults = boost::assign::list_of
        (fLeft1Frp_)
        (fLeft1Ffp_)
        (fRight1Ffp_)
        (fRight1Rfp_)
        (fLeft3FrpMb0_)
        (f1Fsh_)(f3Fsh_);
    testNoDifferences(fInput, fExpectedResults);

    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsInput = boost::assign::list_of
        (rRight1Frp_)(rRight2Frp_)(rRight3Frp_)
        (rLeft1Rrp_)(rLeft2Rrp_)
        (rRight1Rrp_)(rRight2Rrp_)
        (rLeft1Rfp_)(rLeft2Rfp_)(rLeft3Rfp_)
        (rRight1FrpMb1_)(rRight2FrpMb1_)
        (sh1Fsh_)(sh2Fsh_)(sh3Fsh_)(sh4Fsh_);
    std::vector<isaac::build::RStrandOrShadowFragmentIndex> rsExpectedResults = boost::assign::list_of
        (rRight1Frp_)
        (rLeft1Rrp_)
        (rRight1Rrp_)
        (rLeft1Rfp_)
        (sh1Fsh_)(sh3Fsh_);
    testNoDifferences(rsInput, rsExpectedResults);

}
