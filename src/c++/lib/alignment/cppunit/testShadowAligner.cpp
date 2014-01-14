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
#include "testShadowAligner.hh"
#include "BuilderInit.hh"
#include "alignment/TemplateLengthStatistics.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestShadowAligner, registryName("ShadowAligner"));

TestShadowAligner::TestShadowAligner()
    : readMetadataList(getReadMetadataList(81, 92))
    , flowcells(1, isaac::flowcell::Layout("", isaac::flowcell::Layout::Fastq, false, 8, std::vector<unsigned>(),
                                           readMetadataList, isaac::alignment::SeedMetadataList(), "blah"))
    , contigList(getContigList(190, 300, 422))
{
}

void TestShadowAligner::setUp()
{
}

void TestShadowAligner::tearDown()
{
}

static const isaac::alignment::matchSelector::SequencingAdapterList testAdapters;

static const int ELAND_MATCH_SCORE = 2;
static const int ELAND_MISMATCH_SCORE = -1;
static const int ELAND_GAP_OPEN_SCORE = -15;
static const int ELAND_GAP_EXTEND_SCORE = -3;
static const int ELAND_MIN_GAP_EXTEND_SCORE = 25;

void TestShadowAligner::testRescueShadowShortest()
{
    using isaac::alignment::ShadowAligner;
    using isaac::alignment::TemplateLengthStatistics;
    using isaac::alignment::Cluster;
    using isaac::alignment::FragmentMetadata;
    using isaac::alignment::BandedSmithWaterman;

    ShadowAligner shadowAligner(flowcells, 8, false, ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE);
    {
        const TemplateLengthStatistics tls(200, 400, 312, 38, 26, TemplateLengthStatistics::FRp, TemplateLengthStatistics::RFm, -1);
        const std::vector<char> bcl0(getBcl(readMetadataList, contigList, 0, 0, 0, false, true));
        Cluster cluster0(getMaxReadLength(readMetadataList));
        cluster0.init(readMetadataList, bcl0.begin(), 1101, 999, isaac::alignment::ClusterXy(0,0), true, 0);
        FragmentMetadata fragment0, fragment1;
        std::vector<FragmentMetadata> shadowList(50);
        fragment0.cluster = &cluster0;
        fragment0.readIndex = 0;
        fragment0.contigId = 0;
        fragment0.position = 0;
        fragment0.reverse = false;
        // rescue the first read
        CPPUNIT_ASSERT(shadowAligner.rescueShadow(contigList, fragment0, shadowList, readMetadataList, testAdapters, tls, 0));
        fragment1 = shadowList[0];
        CPPUNIT_ASSERT_EQUAL(fragment0.contigId, fragment1.contigId);
        CPPUNIT_ASSERT_EQUAL(fragment0.cluster, fragment1.cluster);
        CPPUNIT_ASSERT_EQUAL((fragment0.readIndex + 1) % 2, fragment1.readIndex);
        CPPUNIT_ASSERT_EQUAL(98L, fragment1.position);
        CPPUNIT_ASSERT_EQUAL(true, fragment1.reverse);
        CPPUNIT_ASSERT_EQUAL(92U, fragment1.observedLength);
        CPPUNIT_ASSERT_EQUAL(0U, fragment1.mismatchCount);
        CPPUNIT_ASSERT_EQUAL(0U, fragment1.cigarOffset);
        CPPUNIT_ASSERT_EQUAL(1U, fragment1.cigarLength);
        CPPUNIT_ASSERT_EQUAL(92U << 4, shadowAligner.getCigarBuffer()[0]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.00920046, fragment1.logProbability, 0.00000001);
        // rescue the second read
        fragment0 = fragment1;
        CPPUNIT_ASSERT(shadowAligner.rescueShadow(contigList, fragment1, shadowList, readMetadataList, testAdapters, tls, 0));
        fragment0 = shadowList[0];
        CPPUNIT_ASSERT_EQUAL(fragment1.contigId, fragment0.contigId);
        CPPUNIT_ASSERT_EQUAL(fragment1.cluster, fragment0.cluster);
        CPPUNIT_ASSERT_EQUAL((fragment1.readIndex + 1) % 2, fragment0.readIndex);
        CPPUNIT_ASSERT_EQUAL(0L, fragment0.position);
        CPPUNIT_ASSERT_EQUAL(false, fragment0.reverse);
        CPPUNIT_ASSERT_EQUAL(81U, fragment0.observedLength);
        CPPUNIT_ASSERT_EQUAL(0U, fragment0.mismatchCount);
        CPPUNIT_ASSERT_EQUAL(0U, fragment0.cigarOffset);
        CPPUNIT_ASSERT_EQUAL(1U, fragment0.cigarLength);
        CPPUNIT_ASSERT_EQUAL(81U << 4, shadowAligner.getCigarBuffer()[0]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.00810041, fragment0.logProbability, 0.00000001);
        // rescue the shadow
        // check the position, mismatches and log probability
    }
    {
        const TemplateLengthStatistics tls(200, 400, 312, 38, 26, TemplateLengthStatistics::FRm, TemplateLengthStatistics::RFp, -1);
        const std::vector<char> bcl0(getBcl(readMetadataList, contigList, 0, 109, 98, true, false));
        Cluster cluster0(getMaxReadLength(readMetadataList));
        cluster0.init(readMetadataList, bcl0.begin(), 1101, 999, isaac::alignment::ClusterXy(0,0), true, 0);
        FragmentMetadata fragment0, fragment1;
        std::vector<FragmentMetadata> shadowList(50);
        fragment0.cluster = &cluster0;
        fragment0.readIndex = 0;
        fragment0.contigId = 0;
        fragment0.position = 0; //109;
        fragment0.reverse = true;
        // rescue the first read
        CPPUNIT_ASSERT(shadowAligner.rescueShadow(contigList, fragment0, shadowList, readMetadataList, testAdapters, tls, 0));
        fragment1 = shadowList[0];
        CPPUNIT_ASSERT_EQUAL(fragment0.contigId, fragment1.contigId);
        CPPUNIT_ASSERT_EQUAL(fragment0.cluster, fragment1.cluster);
        CPPUNIT_ASSERT_EQUAL((fragment0.readIndex + 1) % 2, fragment1.readIndex);
        CPPUNIT_ASSERT_EQUAL(98L, fragment1.position);
        CPPUNIT_ASSERT_EQUAL(false, fragment1.reverse);
        CPPUNIT_ASSERT_EQUAL(92U, fragment1.observedLength);
        CPPUNIT_ASSERT_EQUAL(0U, fragment1.mismatchCount);
        CPPUNIT_ASSERT_EQUAL(0U, fragment1.cigarOffset);
        CPPUNIT_ASSERT_EQUAL(1U, fragment1.cigarLength);
        CPPUNIT_ASSERT_EQUAL(92U << 4, shadowAligner.getCigarBuffer()[0]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.00920046, fragment1.logProbability, 0.00000001);
        // rescue the second read
        fragment0 = fragment1;
        CPPUNIT_ASSERT(shadowAligner.rescueShadow(contigList, fragment1, shadowList, readMetadataList, testAdapters, tls, 0));
        fragment0 = shadowList[0];
        CPPUNIT_ASSERT_EQUAL(fragment1.contigId, fragment0.contigId);
        CPPUNIT_ASSERT_EQUAL(fragment1.cluster, fragment0.cluster);
        CPPUNIT_ASSERT_EQUAL((fragment1.readIndex + 1) % 2, fragment0.readIndex);
        CPPUNIT_ASSERT_EQUAL(0L, fragment0.position);
        CPPUNIT_ASSERT_EQUAL(true, fragment0.reverse);
        CPPUNIT_ASSERT_EQUAL(81U, fragment0.observedLength);
        CPPUNIT_ASSERT_EQUAL(0U, fragment0.mismatchCount);
        CPPUNIT_ASSERT_EQUAL(0U, fragment0.cigarOffset);
        CPPUNIT_ASSERT_EQUAL(1U, fragment0.cigarLength);
        CPPUNIT_ASSERT_EQUAL(81U << 4, shadowAligner.getCigarBuffer()[0]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.00810041, fragment0.logProbability, 0.00000001);
        // rescue the shadow
        // check the position, mismatches and log probability
    }
}

void TestShadowAligner::testRescueShadowLongest()
{
    using isaac::alignment::ShadowAligner;
    using isaac::alignment::TemplateLengthStatistics;
    using isaac::alignment::Cluster;
    using isaac::alignment::FragmentMetadata;
    using isaac::alignment::BandedSmithWaterman;

    ShadowAligner shadowAligner(flowcells, 8, false, ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE);
    {
        const TemplateLengthStatistics tls(200, 400, 312, 38, 26, TemplateLengthStatistics::FRp, TemplateLengthStatistics::RFm, -1);
        const std::vector<char> bcl0(getBcl(readMetadataList, contigList, 4, 0, 12, false, true));
        Cluster cluster0(getMaxReadLength(readMetadataList));
        cluster0.init(readMetadataList, bcl0.begin(), 1101, 999, isaac::alignment::ClusterXy(0,0), true, 0);
        FragmentMetadata fragment0, fragment1;
        std::vector<FragmentMetadata> shadowList(50);
        fragment0.cluster = &cluster0;
        fragment0.readIndex = 0;
        fragment0.contigId = 4;
        fragment0.position = 0;
        fragment0.reverse = false;
        // rescue the first read
        CPPUNIT_ASSERT(shadowAligner.rescueShadow(contigList, fragment0, shadowList, readMetadataList, testAdapters, tls, 0));
        fragment1 = shadowList[0];
        CPPUNIT_ASSERT_EQUAL(fragment0.contigId, fragment1.contigId);
        CPPUNIT_ASSERT_EQUAL(fragment0.cluster, fragment1.cluster);
        CPPUNIT_ASSERT_EQUAL((fragment0.readIndex + 1) % 2, fragment1.readIndex);
        CPPUNIT_ASSERT_EQUAL(318L, fragment1.position);
        CPPUNIT_ASSERT_EQUAL(true, fragment1.reverse);
        CPPUNIT_ASSERT_EQUAL(92U, fragment1.observedLength);
        CPPUNIT_ASSERT_EQUAL(0U, fragment1.mismatchCount);
        CPPUNIT_ASSERT_EQUAL(0U, fragment1.cigarOffset);
        CPPUNIT_ASSERT_EQUAL(1U, fragment1.cigarLength);
        CPPUNIT_ASSERT_EQUAL(92U << 4, shadowAligner.getCigarBuffer()[0]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.00920046, fragment1.logProbability, 0.00000001);
        // rescue the second read
        fragment0 = fragment1;
        CPPUNIT_ASSERT(shadowAligner.rescueShadow(contigList, fragment1, shadowList, readMetadataList, testAdapters, tls, 0));
        fragment0 = shadowList[0];
        CPPUNIT_ASSERT_EQUAL(fragment1.contigId, fragment0.contigId);
        CPPUNIT_ASSERT_EQUAL(fragment1.cluster, fragment0.cluster);
        CPPUNIT_ASSERT_EQUAL((fragment1.readIndex + 1) % 2, fragment0.readIndex);
        CPPUNIT_ASSERT_EQUAL(0L, fragment0.position);
        CPPUNIT_ASSERT_EQUAL(false, fragment0.reverse);
        CPPUNIT_ASSERT_EQUAL(81U, fragment0.observedLength);
        CPPUNIT_ASSERT_EQUAL(0U, fragment0.mismatchCount);
        CPPUNIT_ASSERT_EQUAL(0U, fragment0.cigarOffset);
        CPPUNIT_ASSERT_EQUAL(1U, fragment0.cigarLength);
        CPPUNIT_ASSERT_EQUAL(81U << 4, shadowAligner.getCigarBuffer()[0]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.00810041, fragment0.logProbability, 0.00000001);
        // rescue the shadow
        // check the position, mismatches and log probability
    }
    {
        const TemplateLengthStatistics tls(200, 400, 312, 38, 26, TemplateLengthStatistics::FRm, TemplateLengthStatistics::RFp, -1);
        const std::vector<char> bcl0(getBcl(readMetadataList, contigList, 4, 341, 318, true, false));
        Cluster cluster0(getMaxReadLength(readMetadataList));
        cluster0.init(readMetadataList, bcl0.begin(), 1101, 999, isaac::alignment::ClusterXy(0,0), true, 0);
        FragmentMetadata fragment0, fragment1;
        std::vector<FragmentMetadata> shadowList(50);
        fragment0.cluster = &cluster0;
        fragment0.readIndex = 0;
        fragment0.contigId = 4;
        fragment0.position = 0; //109;
        fragment0.reverse = true;
        // rescue the first read
        CPPUNIT_ASSERT(shadowAligner.rescueShadow(contigList, fragment0, shadowList, readMetadataList, testAdapters, tls, 0));
        fragment1 = shadowList[0];
        CPPUNIT_ASSERT_EQUAL(fragment0.contigId, fragment1.contigId);
        CPPUNIT_ASSERT_EQUAL(fragment0.cluster, fragment1.cluster);
        CPPUNIT_ASSERT_EQUAL((fragment0.readIndex + 1) % 2, fragment1.readIndex);
        CPPUNIT_ASSERT_EQUAL(318L, fragment1.position);
        CPPUNIT_ASSERT_EQUAL(false, fragment1.reverse);
        CPPUNIT_ASSERT_EQUAL(92U, fragment1.observedLength);
        CPPUNIT_ASSERT_EQUAL(0U, fragment1.mismatchCount);
        CPPUNIT_ASSERT_EQUAL(0U, fragment1.cigarOffset);
        CPPUNIT_ASSERT_EQUAL(1U, fragment1.cigarLength);
        CPPUNIT_ASSERT_EQUAL(92U << 4, shadowAligner.getCigarBuffer()[0]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.00920046, fragment1.logProbability, 0.00000001);
        // rescue the second read
        fragment0 = fragment1;
        CPPUNIT_ASSERT(shadowAligner.rescueShadow(contigList, fragment1, shadowList, readMetadataList, testAdapters, tls, 0));
        fragment0 = shadowList[0];
        CPPUNIT_ASSERT_EQUAL(fragment1.contigId, fragment0.contigId);
        CPPUNIT_ASSERT_EQUAL(fragment1.cluster, fragment0.cluster);
        CPPUNIT_ASSERT_EQUAL((fragment1.readIndex + 1) % 2, fragment0.readIndex);
        CPPUNIT_ASSERT_EQUAL(0L, fragment0.position);
        CPPUNIT_ASSERT_EQUAL(true, fragment0.reverse);
        CPPUNIT_ASSERT_EQUAL(81U, fragment0.observedLength);
        CPPUNIT_ASSERT_EQUAL(0U, fragment0.mismatchCount);
        CPPUNIT_ASSERT_EQUAL(0U, fragment0.cigarOffset);
        CPPUNIT_ASSERT_EQUAL(1U, fragment0.cigarLength);
        CPPUNIT_ASSERT_EQUAL(81U << 4, shadowAligner.getCigarBuffer()[0]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.00810041, fragment0.logProbability, 0.00000001);
        // rescue the shadow
        // check the position, mismatches and log probability
    }
}
