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

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>
#include <boost/assign/std/vector.hpp> 

using namespace std;

#include "RegistryName.hh"
#include "testTemplateBuilder.hh"
#include "BuilderInit.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestTemplateBuilder, registryName("TemplateBuilder"));

isaac::alignment::FragmentMetadata getFragmentMetadata(
    unsigned contigId,
    long position,
    long observedLength,
    unsigned readIndex,
    bool reverse,
    unsigned cigarOffset,
    unsigned cigarLength,
    const std::vector<unsigned> *cigarBuffer,
    unsigned mismatchCount,
    double logProbability,
    unsigned uniqueSeedCount,
    unsigned alignmentScore,
    const isaac::alignment::Cluster *cluster)
{
    isaac::alignment::FragmentMetadata f;
    f.contigId = contigId;
    f.position = position;
    f.observedLength = observedLength;
    f.readIndex = readIndex;
    f.reverse = reverse;
    f.cigarOffset = cigarOffset;
    f.cigarLength = cigarLength;
    f.cigarBuffer = cigarBuffer;
    f.mismatchCount = mismatchCount;
    f.logProbability = logProbability;
    f.uniqueSeedCount = uniqueSeedCount;
    f.alignmentScore = alignmentScore;
    f.cluster = cluster;
    return f;
}

TestTemplateBuilder::TestTemplateBuilder()
    : readMetadataList(getReadMetadataList())
    , flowcells(1, isaac::flowcell::Layout("", isaac::flowcell::Layout::Fastq, std::vector<unsigned>(),
                                          readMetadataList, isaac::alignment::SeedMetadataList(), "blah"))
    , contigList(getContigList())
    , cigarBuffer(1000, 1600)
    , tls(readMetadataList, contigList)
    , bcl0(getBcl(readMetadataList, contigList, 0, 2, 3))
    , bcl2(getBcl(readMetadataList, contigList, 2, 1, 2))
    , tile0(32)
    , tile2(31)
    , clusterId0(1234)
    , clusterId2(12345)
    , cluster0(isaac::flowcell::getMaxReadLength(readMetadataList))
    , cluster2(isaac::flowcell::getMaxReadLength(readMetadataList))
      //, f0_0(getFragmentMetadata(0,2,100,0, false, 0, 1, &cigarBuffer, 0, -78.0, 3, 254, &cluster0))
      //, f0_1(getFragmentMetadata(0,107,99,1, true, 1, 1, &cigarBuffer, 2, -92.0, 1, 253, &cluster0))
      , f0_0(getFragmentMetadata(0,2,100,0, false, 0, 1, &cigarBuffer, 0, -8.0, 3, 254, &cluster0))
      , f0_1(getFragmentMetadata(0,107,99,1, true, 1, 1, &cigarBuffer, 2, -12.0, 1, 253, &cluster0))
{
    cluster0.init(readMetadataList, bcl0.begin(), tile0, clusterId0, isaac::alignment::ClusterXy(0,0), true, 0);
    cluster2.init(readMetadataList, bcl2.begin(), tile2, clusterId2, isaac::alignment::ClusterXy(0,0), true, 0);
}

void TestTemplateBuilder::setUp()
{
}

void TestTemplateBuilder::tearDown()
{
}

void TestTemplateBuilder::testConstructor()
{
    CPPUNIT_ASSERT_EQUAL(150U, tls.getMin());
    CPPUNIT_ASSERT_EQUAL(190U, tls.getMedian());
    CPPUNIT_ASSERT_EQUAL(250U, tls.getMax());
    CPPUNIT_ASSERT_EQUAL(20U, tls.getLowStdDev());
    CPPUNIT_ASSERT_EQUAL(30U, tls.getHighStdDev());
    CPPUNIT_ASSERT(tls.isStable());
    CPPUNIT_ASSERT_EQUAL(std::string("FR+"), tls.alignmentModelName(tls.getBestModel(0)));
    CPPUNIT_ASSERT_EQUAL(std::string("RF-"), tls.alignmentModelName(tls.getBestModel(1)));
    CPPUNIT_ASSERT_EQUAL(1U, readMetadataList[1].getIndex());
    CPPUNIT_ASSERT_EQUAL(1U, cluster0[1].getIndex());
    CPPUNIT_ASSERT_EQUAL(254U, f0_0.alignmentScore);
}

void TestTemplateBuilder::checkUnalignedTemplate(
    const isaac::alignment::BamTemplate &bamTemplate,
    const isaac::alignment::Cluster &cluster) const
{
    checkUnalignedFragment(bamTemplate, cluster, 0, 0);
    checkUnalignedFragment(bamTemplate, cluster, 1, 1);
}


void TestTemplateBuilder::checkUnalignedFragment(
    const isaac::alignment::BamTemplate &bamTemplate,
    const isaac::alignment::Cluster &cluster,
    const unsigned i,
    const unsigned readIndex) const
{
    CPPUNIT_ASSERT_EQUAL(i, readIndex);
    // check for an unaligned fragment
    CPPUNIT_ASSERT(bamTemplate.getFragmentMetadata(i).isNoMatch());
    CPPUNIT_ASSERT_EQUAL(0L, bamTemplate.getFragmentMetadata(i).observedLength);
    CPPUNIT_ASSERT_EQUAL(readIndex, bamTemplate.getFragmentMetadata(i).readIndex);
    CPPUNIT_ASSERT_EQUAL(false, bamTemplate.getFragmentMetadata(i).reverse);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(i).cigarOffset);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(i).cigarLength);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(i).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(0.0, bamTemplate.getFragmentMetadata(i).logProbability);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(i).uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL(-1U, bamTemplate.getFragmentMetadata(i).alignmentScore);
    CPPUNIT_ASSERT_EQUAL(&cluster, bamTemplate.getFragmentMetadata(i).cluster);
}

static const isaac::alignment::matchSelector::SequencingAdapterList testAdapters;

static const int ELAND_MATCH_SCORE = 2;
static const int ELAND_MISMATCH_SCORE = -1;
static const int ELAND_GAP_OPEN_SCORE = -15;
static const int ELAND_GAP_EXTEND_SCORE = -3;
static const int ELAND_MIN_GAP_EXTEND_SCORE = 25;

void TestTemplateBuilder::testEmptyMatchList()
{
    using isaac::alignment::TemplateBuilder;
    using isaac::alignment::BamTemplate;
    using isaac::alignment::FragmentMetadata;
    using isaac::alignment::BandedSmithWaterman;
    std::auto_ptr<TemplateBuilder> templateBuilder(new TemplateBuilder(flowcells, 10, 4, false, 8, false,
                                                                       ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE,
                                                                       ELAND_MIN_GAP_EXTEND_SCORE, 20000,
                                                                       TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED));
    const BamTemplate &bamTemplate = templateBuilder->getBamTemplate();
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentCount());
    std::vector<std::vector<FragmentMetadata> > fragments(2);
    templateBuilder->buildTemplate(contigList, readMetadataList, testAdapters, fragments, cluster0, tls);
    CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getFragmentCount());
    checkUnalignedTemplate(bamTemplate, cluster0);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getAlignmentScore());
    // initialize the first fragment with garbage
    fragments[0].push_back(f0_0);
    fragments[1].push_back(f0_1);
    templateBuilder->buildTemplate(contigList, readMetadataList, testAdapters, fragments, cluster0, tls);
    // clear the fragments and check again
    fragments[0].clear();
    fragments[1].clear();
    templateBuilder->buildTemplate(contigList, readMetadataList, testAdapters, fragments, cluster0, tls);
    CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getFragmentCount());
    checkUnalignedTemplate(bamTemplate, cluster0);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getAlignmentScore());
}

void TestTemplateBuilder::testOrphan()
{
    using isaac::alignment::TemplateBuilder;
    using isaac::alignment::BamTemplate;
    using isaac::alignment::FragmentMetadata;
    using isaac::alignment::BandedSmithWaterman;
    std::auto_ptr<TemplateBuilder> templateBuilder(new TemplateBuilder(flowcells, 10, 4, false, 8, false,
                                                                       ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE,
                                                                       ELAND_MIN_GAP_EXTEND_SCORE, 20000,
                                                                       TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED));
    const BamTemplate &bamTemplate = templateBuilder->getBamTemplate();
    std::vector<std::vector<FragmentMetadata> > fragments(2);
    // align on the first read only
    fragments[0].push_back(f0_0);
    templateBuilder->buildTemplate(contigList, readMetadataList, testAdapters, fragments, cluster0, tls);
    // this orphan should be rescued
    CPPUNIT_ASSERT_EQUAL(1136U, bamTemplate.getAlignmentScore());
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).contigId);
    CPPUNIT_ASSERT_EQUAL(2L, bamTemplate.getFragmentMetadata(0).position);
    CPPUNIT_ASSERT_EQUAL(100L, bamTemplate.getFragmentMetadata(0).observedLength);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).readIndex);
    CPPUNIT_ASSERT_EQUAL(false, bamTemplate.getFragmentMetadata(0).reverse);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).cigarOffset);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(0).cigarLength);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(f0_0.logProbability, bamTemplate.getFragmentMetadata(0).logProbability);
    CPPUNIT_ASSERT_EQUAL(3U, bamTemplate.getFragmentMetadata(0).uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL(534U, bamTemplate.getFragmentMetadata(0).alignmentScore);
    CPPUNIT_ASSERT_EQUAL(569U, bamTemplate.getFragmentMetadata(1).alignmentScore);
    CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(0).cluster);
    // align on the second read only
    fragments[0].clear();
    fragments[1].push_back(f0_1);
    templateBuilder->buildTemplate(contigList, readMetadataList, testAdapters, fragments, cluster0, tls);
    // this one should be rescued as well
    CPPUNIT_ASSERT_EQUAL(1119U, bamTemplate.getAlignmentScore());
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(1).contigId);
    CPPUNIT_ASSERT_EQUAL(107L, bamTemplate.getFragmentMetadata(1).position);
    CPPUNIT_ASSERT_EQUAL(99L, bamTemplate.getFragmentMetadata(1).observedLength);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).readIndex);
    CPPUNIT_ASSERT_EQUAL(true, bamTemplate.getFragmentMetadata(1).reverse);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).cigarOffset);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).cigarLength);
    CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getFragmentMetadata(1).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(f0_1.logProbability, bamTemplate.getFragmentMetadata(1).logProbability);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL(517U, bamTemplate.getFragmentMetadata(1).alignmentScore);
    CPPUNIT_ASSERT_EQUAL(569U, bamTemplate.getFragmentMetadata(0).alignmentScore);
    CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(1).cluster);
}

void TestTemplateBuilder::testUnique()
{
    using isaac::alignment::TemplateBuilder;
    using isaac::alignment::BamTemplate;
    using isaac::alignment::FragmentMetadata;
    using isaac::alignment::BandedSmithWaterman;
    std::auto_ptr<TemplateBuilder> templateBuilder(new TemplateBuilder(flowcells, 10, 4, false, 8, false,
                                                                       ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE,
                                                                       ELAND_MIN_GAP_EXTEND_SCORE, 20000,
                                                                       TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED));
    const BamTemplate &bamTemplate = templateBuilder->getBamTemplate();
    std::vector<std::vector<FragmentMetadata> > fragments(2);
    fragments[0].push_back(f0_0);
    fragments[1].push_back(f0_1);
    templateBuilder->buildTemplate(contigList, readMetadataList, testAdapters, fragments, cluster0, tls);
    CPPUNIT_ASSERT_EQUAL(1084U, bamTemplate.getAlignmentScore());
    //CPPUNIT_ASSERT_EQUAL(bamTemplate.getFragmentMetadata(0).getAlignmentScore() + bamTemplate.getFragmentMetadata(1).getAlignmentScore(),
    //                     bamTemplate.getAlignmentScore());
    CPPUNIT_ASSERT_EQUAL(534U, bamTemplate.getFragmentMetadata(0).getAlignmentScore());
    CPPUNIT_ASSERT_EQUAL(517U, bamTemplate.getFragmentMetadata(1).getAlignmentScore());
    // check the first read
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).contigId);
    CPPUNIT_ASSERT_EQUAL(2L, bamTemplate.getFragmentMetadata(0).position);
    CPPUNIT_ASSERT_EQUAL(100L, bamTemplate.getFragmentMetadata(0).observedLength);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).readIndex);
    CPPUNIT_ASSERT_EQUAL(false, bamTemplate.getFragmentMetadata(0).reverse);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).cigarOffset);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(0).cigarLength);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(f0_0.logProbability, bamTemplate.getFragmentMetadata(0).logProbability);
    CPPUNIT_ASSERT_EQUAL(3U, bamTemplate.getFragmentMetadata(0).uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL(534U, bamTemplate.getFragmentMetadata(0).alignmentScore);
    CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(0).cluster);
    // check the second read
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(1).contigId);
    CPPUNIT_ASSERT_EQUAL(107L, bamTemplate.getFragmentMetadata(1).position);
    CPPUNIT_ASSERT_EQUAL(99L, bamTemplate.getFragmentMetadata(1).observedLength);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).readIndex);
    CPPUNIT_ASSERT_EQUAL(true, bamTemplate.getFragmentMetadata(1).reverse);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).cigarOffset);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).cigarLength);
    CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getFragmentMetadata(1).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(f0_1.logProbability, bamTemplate.getFragmentMetadata(1).logProbability);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL(517U, bamTemplate.getFragmentMetadata(1).alignmentScore);
    CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(1).cluster);
}

/**
 * \brief This test was originally designed to ensure the pair that matches the tls is picked.
 * Currently, everything on the same contig with the correct orientation with size below max_ + 50000
 * is considered allright. So, verify that the best alignment score pair is picked.
 */
void TestTemplateBuilder::testMultiple()
{
    using isaac::alignment::TemplateBuilder;
    using isaac::alignment::BamTemplate;
    using isaac::alignment::FragmentMetadata;
    using isaac::alignment::BandedSmithWaterman;
    std::auto_ptr<TemplateBuilder> templateBuilder(new TemplateBuilder(flowcells, 10, 4, false, 8, false,
                                                                       ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE,
                                                                       ELAND_MIN_GAP_EXTEND_SCORE, 20000,
                                                                       TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED));
    const BamTemplate &bamTemplate = templateBuilder->getBamTemplate();

    std::vector<std::vector<FragmentMetadata> > fragments(2);
    FragmentMetadata t0 = f0_0;
    FragmentMetadata t1 = f0_1;
    for (size_t i = 0; 2 > i; ++i)
    {
        fragments[0].push_back(t0);
        t0.position += 56;
        fragments[0].push_back(t0);
        t0.position += 65;
        fragments[1].push_back(t1);
        t1.position += 300;
    }
    t0 = f0_0;
    t1 = f0_1;
    t0.contigId = 1;
    t1.contigId = 1;
    for (size_t i = 0; 2 > i; ++i)
    {
        t0.position += 56;
        fragments[0].push_back(t0);
        t0.position += 65;
        fragments[0].push_back(t0);
        t1.position += 401;
        fragments[1].push_back(t1);
    }
    t0 = f0_0;
    t1 = f0_1;
    t0.contigId = 1;
    t1.contigId = 1;
    t0.logProbability += 2;
    t1.logProbability += 2;
    fragments[0].push_back(t0);
    FragmentMetadata best0 = fragments[0].back();
    fragments[1].push_back(t1);
    FragmentMetadata best1 = fragments[1].back();
    t0.logProbability -= 2;
    t1.logProbability -= 2;
    for (size_t i = 0; 2 > i; ++i)
    {
        t0.position += 36;
        fragments[0].push_back(t0);
        t0.position += 45;
        fragments[0].push_back(t0);
        t1.position += 402;
        fragments[1].push_back(t1);
    }
    templateBuilder->buildTemplate(contigList, readMetadataList, testAdapters, fragments, cluster0, tls);

    CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getAlignmentScore());
    // check the first read
    CPPUNIT_ASSERT_EQUAL(unsigned(best0.getFStrandReferencePosition().getContigId()), bamTemplate.getFragmentMetadata(0).contigId);
    CPPUNIT_ASSERT_EQUAL(long(best0.getFStrandReferencePosition().getPosition()), bamTemplate.getFragmentMetadata(0).position);
    CPPUNIT_ASSERT_EQUAL(best0.getObservedLength(), unsigned(bamTemplate.getFragmentMetadata(0).observedLength));
    CPPUNIT_ASSERT_EQUAL(best0.getReadIndex(), bamTemplate.getFragmentMetadata(0).readIndex);
    CPPUNIT_ASSERT_EQUAL(best0.isReverse(), bamTemplate.getFragmentMetadata(0).reverse);
    CPPUNIT_ASSERT_EQUAL(best0.cigarOffset, bamTemplate.getFragmentMetadata(0).cigarOffset);
    CPPUNIT_ASSERT_EQUAL(best0.getCigarLength(), bamTemplate.getFragmentMetadata(0).cigarLength);
    CPPUNIT_ASSERT_EQUAL(best0.getMismatchCount(), bamTemplate.getFragmentMetadata(0).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(best0.logProbability, bamTemplate.getFragmentMetadata(0).logProbability);
    CPPUNIT_ASSERT_EQUAL(best0.uniqueSeedCount, bamTemplate.getFragmentMetadata(0).uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getFragmentMetadata(0).alignmentScore);
    CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(0).cluster);
    // check the second read
    CPPUNIT_ASSERT_EQUAL(unsigned(best1.getFStrandReferencePosition().getContigId()), bamTemplate.getFragmentMetadata(1).contigId);
    CPPUNIT_ASSERT_EQUAL(long(best1.getFStrandReferencePosition().getPosition()), bamTemplate.getFragmentMetadata(1).position);
    CPPUNIT_ASSERT_EQUAL(best1.getObservedLength(), unsigned(bamTemplate.getFragmentMetadata(1).observedLength));
    CPPUNIT_ASSERT_EQUAL(best1.getReadIndex(), bamTemplate.getFragmentMetadata(1).readIndex);
    CPPUNIT_ASSERT_EQUAL(best1.isReverse(), bamTemplate.getFragmentMetadata(1).reverse);
    CPPUNIT_ASSERT_EQUAL(best1.cigarOffset, bamTemplate.getFragmentMetadata(1).cigarOffset);
    CPPUNIT_ASSERT_EQUAL(best1.getCigarLength(), bamTemplate.getFragmentMetadata(1).cigarLength);
    CPPUNIT_ASSERT_EQUAL(best1.getMismatchCount(), bamTemplate.getFragmentMetadata(1).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(best1.logProbability, bamTemplate.getFragmentMetadata(1).logProbability);
    CPPUNIT_ASSERT_EQUAL(best1.uniqueSeedCount, bamTemplate.getFragmentMetadata(1).uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL(3U, bamTemplate.getFragmentMetadata(1).alignmentScore);
    CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(1).cluster);

}

DummyTemplateLengthStatistics::DummyTemplateLengthStatistics(
    const std::vector<isaac::flowcell::ReadMetadata> readMetadataList,
    const std::vector<isaac::reference::Contig> contigList) :
    TemplateLengthStatistics(-1)
{
    reset(contigList, readMetadataList);

    setMin(150);
    setMax(250);
    setMedian(190);
    setLowStdDev(20);
    setHighStdDev(30);
    setBestModel(FRp, 0); // FR+
    setBestModel(RFm, 1); // RF-
    setStable(true);
}
