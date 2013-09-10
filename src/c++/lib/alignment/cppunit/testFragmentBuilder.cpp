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
#include "testFragmentBuilder.hh"
#include "BuilderInit.hh"

#include "flowcell/SequencingAdapterMetadata.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestFragmentBuilder, registryName("FragmentBuilder"));


TestFragmentBuilder::TestFragmentBuilder()
    : readMetadataList(getReadMetadataList())
    , seedMetadataList(getSeedMetadataList())
    , flowcells(1, isaac::flowcell::Layout("", isaac::flowcell::Layout::Fastq, std::vector<unsigned>(),
                                           readMetadataList, seedMetadataList, "blah"))

    , contigList(getContigList())
    , bcl0(getBcl(readMetadataList, contigList, 0, 2, 3))
    , bcl2(getBcl(readMetadataList, contigList, 2, 1, 2))
    , bcl3(subv(bcl0, 0,4) +
           "x" +  // x == 30*4+0 -- replaces 'T' with 'A'
           subv(bcl0, 5, 190) +
           "Q" +  // x == 20*4+0 -- replaces 'A' with 'C'
           subv(bcl0, 196))
    , bcl4l(getBcl(substr(contigList[4].forward_, 0, 44) + substr(contigList[4].forward_, 0, 56) +
                   substr(reverseComplement(contigList[4].forward_), 0, 42) + substr(reverseComplement(contigList[4].forward_), 0, 58)))
    , bcl4t(getBcl(substr(contigList[4].forward_, 16, 44) + substr(contigList[4].forward_, 0, 56) +
                   substr(reverseComplement(contigList[4].forward_), 18, 42) + substr(reverseComplement(contigList[4].forward_), 0, 58)))
    , bcl4lt(getBcl(std::vector<char>(10, 'A') + contigList[4].forward_ + std::vector<char>(30, 'C') +
                    std::vector<char>(15, 'G') + reverseComplement(contigList[4].forward_) + std::vector<char>(25, 'T')))
    , tile0(32)
    , tile2(31)
    , clusterId0(1234)
    , clusterId2(12345)
    , cluster0(getMaxReadLength(readMetadataList))
    , cluster2(getMaxReadLength(readMetadataList))
    , cluster3(getMaxReadLength(readMetadataList))
    , cluster4l(getMaxReadLength(readMetadataList))
    , cluster4t(getMaxReadLength(readMetadataList))
    , cluster4lt(getMaxReadLength(readMetadataList))
{
    cluster0.init(readMetadataList, bcl0.begin(), tile0, clusterId0, isaac::alignment::ClusterXy(0,0), true, 0);
    cluster2.init(readMetadataList, bcl2.begin(), tile2, clusterId2, isaac::alignment::ClusterXy(0,0), true, 0);
    cluster3.init(readMetadataList, bcl3.begin(), tile0, clusterId0, isaac::alignment::ClusterXy(0,0), true, 0);
    cluster4l.init(readMetadataList, bcl4l.begin(), tile0, clusterId0, isaac::alignment::ClusterXy(0,0), true, 0);
    cluster4t.init(readMetadataList, bcl4t.begin(), tile0, clusterId0, isaac::alignment::ClusterXy(0,0), true, 0);
    cluster4lt.init(readMetadataList, bcl4lt.begin(), tile0, clusterId0, isaac::alignment::ClusterXy(0,0), true, 0);
}

void TestFragmentBuilder::setUp()
{
    matchList.clear();
}

void TestFragmentBuilder::tearDown()
{
}

static const isaac::alignment::matchSelector::SequencingAdapterList testAdapters;

static const int ELAND_MATCH_SCORE = 2;
static const int ELAND_MISMATCH_SCORE = -1;
static const int ELAND_GAP_OPEN_SCORE = -15;
static const int ELAND_GAP_EXTEND_SCORE = -3;
static const int ELAND_MIN_GAP_EXTEND_SCORE = 25;

void TestFragmentBuilder::testEmptyMatchList()
{
    using isaac::alignment::FragmentBuilder;
    FragmentBuilder fragmentBuilder(flowcells, 123, seedMetadataList.size()/2, 8, false,
                                    ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE,
                                    ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    // check the emptyness after creation
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragmentBuilder.getFragments().size());
    CPPUNIT_ASSERT(fragmentBuilder.getFragments()[0].empty());
    CPPUNIT_ASSERT(fragmentBuilder.getFragments()[1].empty());
    CPPUNIT_ASSERT(fragmentBuilder.getCigarBuffer().empty());
    // build fragments for an empty list
    fragmentBuilder.build(contigList, readMetadataList, seedMetadataList, testAdapters,
                          matchList.begin(), matchList.end(), cluster0, true);
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragmentBuilder.getFragments().size());
    CPPUNIT_ASSERT(fragmentBuilder.getFragments()[0].empty());
    CPPUNIT_ASSERT(fragmentBuilder.getFragments()[1].empty());
    CPPUNIT_ASSERT(fragmentBuilder.getCigarBuffer().empty());
}

void TestFragmentBuilder::auxSingleSeed(const unsigned s0, const unsigned s1)
{
    const unsigned offset0 = seedMetadataList[s0].getOffset();
    const unsigned offset1 = seedMetadataList[s1].getOffset();
    // build fragments for a single seed on each read
    // SeedId(tile, cluster, seed, reverse)
    // ReferencePosition(contigId, position)
    using isaac::alignment::SeedId;
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, s0, false), ReferencePosition(0, 2 + offset0))));
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, s1, true ), ReferencePosition(0, 175 - offset1))));
    // Sanity check on the match list
    CPPUNIT_ASSERT_EQUAL((size_t)2, matchList.size());
    CPPUNIT_ASSERT_EQUAL((unsigned long)s0, matchList[0].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)s1, matchList[1].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, seedMetadataList[matchList[0].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)1, seedMetadataList[matchList[1].seedId.getSeed()].getReadIndex());
    // Create the fragment builder
    using isaac::alignment::FragmentBuilder;
    using isaac::alignment::Cigar;
    FragmentBuilder fragmentBuilder(flowcells, 456, seedMetadataList.size()/2, 8, false,
                                    ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE,
                                    ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    // build the fragments
    fragmentBuilder.build(contigList, readMetadataList, seedMetadataList, testAdapters, matchList.begin(), matchList.end(), cluster0, true);
    // check buffer geometry
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragmentBuilder.getFragments().size());
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragmentBuilder.getCigarBuffer().size());
    // Fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((size_t)1, fragmentBuilder.getFragments()[0].size());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].contigId);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[0][0].uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL((long)2, fragmentBuilder.getFragments()[0][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragmentBuilder.getFragments()[0][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragmentBuilder.getFragments()[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragmentBuilder.getFragments()[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragmentBuilder.getFragments()[0][0].logProbability, (double)0.000001);
    // Check that the fragments actually align at the expected positions
    for (size_t i = 0; 100 > i; ++i)
    {
        const unsigned position = fragmentBuilder.getFragments()[0][0].position;
        CPPUNIT_ASSERT_EQUAL(cluster0[0].getForwardSequence()[i], contigList[0].forward_[i + position]);
    }
    // Fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((size_t)1, fragmentBuilder.getFragments()[1].size());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[1][0].contigId);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[1][0].uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL((long)107, fragmentBuilder.getFragments()[1][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragmentBuilder.getFragments()[1][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragmentBuilder.getFragments()[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragmentBuilder.getFragments()[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[1][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragmentBuilder.getFragments()[1][0].logProbability, (double)0.000001);
    // Check that the fragments actually align at the expected positions
    const std::vector<char> reverse0 = reverseComplement(contigList[0].forward_);
    for (size_t i = 0; 100 > i; ++i)
    {
        CPPUNIT_ASSERT_EQUAL(cluster0[1].getForwardSequence()[i], reverse0[i + 3]);
        const unsigned position = fragmentBuilder.getFragments()[1][0].position;
        const unsigned length = fragmentBuilder.getFragments()[1][0].observedLength;
        const char reverseBase = isaac::oligo::getReverseBase(cluster0[1].getForwardSequence()[i]);
        CPPUNIT_ASSERT_EQUAL(reverseBase, contigList[0].forward_[position + length - 1 - i]);
    }
}

void TestFragmentBuilder::testSingleSeed()
{
    // test on seed index 0 and 3 (both at offset 0)
    auxSingleSeed(0, 3);
}

void TestFragmentBuilder::testSeedOffset()
{
    // test on seed index 0 and 3 (both at offset 0)
    auxSingleSeed(1, 5);
}

void TestFragmentBuilder::testMultiSeed()
{
    // build fragments for a single seed on each read
    // SeedId(tile, cluster, seed, reverse)
    // ReferencePosition(contigId, position)
    using isaac::alignment::SeedId;
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, 0, false), ReferencePosition(0, 2))));
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, 1, false), ReferencePosition(0, 2 + 32))));
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, 2, false), ReferencePosition(0, 2 + 64))));
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, 3, true ), ReferencePosition(0, 175))));
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, 4, true ), ReferencePosition(0, 175 - 32))));
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, 5, true ), ReferencePosition(0, 175 - 64))));
    // Sanity check on the match list
    CPPUNIT_ASSERT_EQUAL((size_t)6, matchList.size());
    CPPUNIT_ASSERT_EQUAL((unsigned long)0, matchList[0].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)1, matchList[1].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)2, matchList[2].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)3, matchList[3].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)4, matchList[4].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)5, matchList[5].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, seedMetadataList[matchList[0].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, seedMetadataList[matchList[1].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, seedMetadataList[matchList[2].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)1, seedMetadataList[matchList[3].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)1, seedMetadataList[matchList[4].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)1, seedMetadataList[matchList[5].seedId.getSeed()].getReadIndex());
    // Create the fragment builder
    using isaac::alignment::FragmentBuilder;
    using isaac::alignment::Cigar;
    FragmentBuilder fragmentBuilder(flowcells, 123, seedMetadataList.size()/2, 8, false,
                                    ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE,
                                    ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    // build the fragments
    fragmentBuilder.build(contigList, readMetadataList, seedMetadataList, testAdapters, matchList.begin(), matchList.end(), cluster0, true);
    // check buffer geometry
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragmentBuilder.getFragments().size());
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragmentBuilder.getCigarBuffer().size());
    // Fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((size_t)1, fragmentBuilder.getFragments()[0].size());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].contigId);
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragmentBuilder.getFragments()[0][0].uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL((long)2, fragmentBuilder.getFragments()[0][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragmentBuilder.getFragments()[0][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragmentBuilder.getFragments()[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragmentBuilder.getFragments()[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragmentBuilder.getFragments()[0][0].logProbability, (double)0.000001);
    // Fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((size_t)1, fragmentBuilder.getFragments()[1].size());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[1][0].contigId);
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragmentBuilder.getFragments()[1][0].uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL((long)107, fragmentBuilder.getFragments()[1][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragmentBuilder.getFragments()[1][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragmentBuilder.getFragments()[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragmentBuilder.getFragments()[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[1][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragmentBuilder.getFragments()[1][0].logProbability, (double)0.000001);
}

void TestFragmentBuilder::testRepeats()
{
    // build fragments for a single seed on each read
    // SeedId(tile, cluster, seed, reverse)
    // ReferencePosition(contigId, position)
    using isaac::alignment::SeedId;
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(Match(SeedId(tile2, 0, clusterId2, 0, false), ReferencePosition(2, 1))));
    matchList.push_back(Match(Match(SeedId(tile2, 0, clusterId2, 1, false), ReferencePosition(2, 1 + 32))));
    matchList.push_back(Match(Match(SeedId(tile2, 0, clusterId2, 2, false), ReferencePosition(2, 1 + 64))));
    matchList.push_back(Match(Match(SeedId(tile2, 0, clusterId2, 3, true ), ReferencePosition(2, 196))));
    matchList.push_back(Match(Match(SeedId(tile2, 0, clusterId2, 4, true ), ReferencePosition(2, 196 - 32))));
    matchList.push_back(Match(Match(SeedId(tile2, 0, clusterId2, 5, true ), ReferencePosition(2, 196 - 64))));
    matchList.push_back(Match(Match(SeedId(tile2, 0, clusterId2, 0, false), ReferencePosition(3, 6))));
    matchList.push_back(Match(Match(SeedId(tile2, 0, clusterId2, 2, false), ReferencePosition(3, 6 + 64))));
    matchList.push_back(Match(Match(SeedId(tile2, 0, clusterId2, 3, true ), ReferencePosition(3, 201))));
    matchList.push_back(Match(Match(SeedId(tile2, 0, clusterId2, 4, true ), ReferencePosition(3, 201 - 32))));
    // Sanity check on the match list
    CPPUNIT_ASSERT_EQUAL((size_t)10, matchList.size());
    CPPUNIT_ASSERT_EQUAL((unsigned long)0, matchList[0].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)1, matchList[1].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)2, matchList[2].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)3, matchList[3].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)4, matchList[4].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)5, matchList[5].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)0, matchList[6].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)2, matchList[7].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)3, matchList[8].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)4, matchList[9].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, seedMetadataList[matchList[0].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, seedMetadataList[matchList[1].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, seedMetadataList[matchList[2].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)1, seedMetadataList[matchList[3].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)1, seedMetadataList[matchList[4].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)1, seedMetadataList[matchList[5].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, seedMetadataList[matchList[6].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, seedMetadataList[matchList[7].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)1, seedMetadataList[matchList[8].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)1, seedMetadataList[matchList[9].seedId.getSeed()].getReadIndex());
    // Create the fragment builder
    using isaac::alignment::FragmentBuilder;
    using isaac::alignment::Cigar;
    FragmentBuilder fragmentBuilder(flowcells, 123, seedMetadataList.size()/2, 8, false,
                                    ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE,
                                    ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    // build the fragments
    fragmentBuilder.build(contigList, readMetadataList, seedMetadataList, testAdapters, matchList.begin(), matchList.end(), cluster2, true);
    // check buffer geometry
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragmentBuilder.getFragments().size());
    CPPUNIT_ASSERT_EQUAL((size_t)4, fragmentBuilder.getCigarBuffer().size());
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragmentBuilder.getFragments()[0].size());
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragmentBuilder.getFragments()[1].size());
    // First fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragmentBuilder.getFragments()[0][0].contigId);
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragmentBuilder.getFragments()[0][0].uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL((long)1, fragmentBuilder.getFragments()[0][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragmentBuilder.getFragments()[0][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragmentBuilder.getFragments()[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragmentBuilder.getFragments()[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragmentBuilder.getFragments()[0][0].logProbability, (double)0.000001);
    // Second fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragmentBuilder.getFragments()[0][1].contigId);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragmentBuilder.getFragments()[0][1].uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL((long)6, fragmentBuilder.getFragments()[0][1].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragmentBuilder.getFragments()[0][1].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragmentBuilder.getFragments()[0][1].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragmentBuilder.getFragments()[0][1].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[0][1].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[0][1].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][1].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragmentBuilder.getFragments()[0][1].logProbability, (double)0.000001);
    // First fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragmentBuilder.getFragments()[1][0].contigId);
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragmentBuilder.getFragments()[1][0].uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL((long)128, fragmentBuilder.getFragments()[1][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragmentBuilder.getFragments()[1][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragmentBuilder.getFragments()[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragmentBuilder.getFragments()[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragmentBuilder.getFragments()[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[1][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragmentBuilder.getFragments()[1][0].logProbability, (double)0.000001);
    // Seconf fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragmentBuilder.getFragments()[1][1].contigId);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragmentBuilder.getFragments()[1][1].uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL((long)133, fragmentBuilder.getFragments()[1][1].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragmentBuilder.getFragments()[1][1].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragmentBuilder.getFragments()[1][1].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragmentBuilder.getFragments()[1][1].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragmentBuilder.getFragments()[1][1].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[1][1].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[1][1].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragmentBuilder.getFragments()[1][1].logProbability, (double)0.000001);
}

void TestFragmentBuilder::testMismatches()
{
    CPPUNIT_ASSERT_EQUAL(bcl0.size(), bcl3.size());
    // Check the mismatch on the forward strand
    CPPUNIT_ASSERT_EQUAL(bcl0[3], bcl3[3]);
    CPPUNIT_ASSERT((bcl0[4]&3) != (bcl3[4]&3));
    CPPUNIT_ASSERT_EQUAL(bcl0[5], bcl3[5]);
    // Check the mismatch on the reverse strand
    CPPUNIT_ASSERT_EQUAL(bcl0[194], bcl3[194]);
    CPPUNIT_ASSERT((bcl0[195]&3) != (bcl3[195]&3));
    CPPUNIT_ASSERT_EQUAL(bcl0[196], bcl3[196]);
    const unsigned s0 = 0;
    const unsigned s1 = 3;
    const unsigned offset0 = seedMetadataList[s0].getOffset();
    const unsigned offset1 = seedMetadataList[s1].getOffset();
    // build fragments for a single seed on each read
    // SeedId(tile, cluster, seed, reverse)
    // ReferencePosition(contigId, position)
    using isaac::alignment::SeedId;
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, s0, false), ReferencePosition(0, 2 + offset0))));
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, s1, true ), ReferencePosition(0, 175 - offset1))));
    // Sanity check on the match list
    CPPUNIT_ASSERT_EQUAL((size_t)2, matchList.size());
    CPPUNIT_ASSERT_EQUAL((unsigned long)s0, matchList[0].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned long)s1, matchList[1].seedId.getSeed());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, seedMetadataList[matchList[0].seedId.getSeed()].getReadIndex());
    CPPUNIT_ASSERT_EQUAL((unsigned)1, seedMetadataList[matchList[1].seedId.getSeed()].getReadIndex());
    // Create the fragment builder
    using isaac::alignment::FragmentBuilder;
    using isaac::alignment::Cigar;
    FragmentBuilder fragmentBuilder(flowcells, 123, seedMetadataList.size()/2, 8, false,
                                    ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE,
                                    ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    // build the fragments
    fragmentBuilder.build(contigList, readMetadataList, seedMetadataList, testAdapters, matchList.begin(), matchList.end(), cluster3, true);
    // check buffer geometry
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragmentBuilder.getFragments().size());
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragmentBuilder.getCigarBuffer().size());
    // Fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((size_t)1, fragmentBuilder.getFragments()[0].size());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].contigId);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[0][0].uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL((long)2, fragmentBuilder.getFragments()[0][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragmentBuilder.getFragments()[0][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragmentBuilder.getFragments()[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragmentBuilder.getFragments()[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[0][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-8.016268063, fragmentBuilder.getFragments()[0][0].logProbability, (double)0.000000001);
    // Check that the fragments actually align at the expected positions
    CPPUNIT_ASSERT_EQUAL(cluster3[0].getForwardSequence()[4], 'A');
    for (size_t i = 0; 100 > i; ++i)
    {
        if (4 != i)
        {
            CPPUNIT_ASSERT_EQUAL(cluster3[0].getForwardSequence()[i], contigList[0].forward_[i + 2]);
        }
    }
    // Fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((size_t)1, fragmentBuilder.getFragments()[1].size());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[1][0].contigId);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[1][0].uniqueSeedCount);
    CPPUNIT_ASSERT_EQUAL((long)107, fragmentBuilder.getFragments()[1][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragmentBuilder.getFragments()[1][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragmentBuilder.getFragments()[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragmentBuilder.getFragments()[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragmentBuilder.getFragments()[1][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-5.713682970, fragmentBuilder.getFragments()[1][0].logProbability, (double)0.000000001);
    // Check that the fragments actually align at the expected positions
    const std::vector<char> reverse0 = reverseComplement(contigList[0].forward_);
    CPPUNIT_ASSERT(cluster3[1].getForwardSequence()[95] != reverse0[95 + 3]);
    for (size_t i = 0; 100 > i; ++i)
    {
        if (95 != i)
        {
            CPPUNIT_ASSERT_EQUAL(cluster3[1].getForwardSequence()[i], reverse0[i + 3]);
        }
    }
}

void TestFragmentBuilder::testLeadingSoftClips()
{
    const unsigned s0 = 2;
    const unsigned s1 = 5;
    // build fragments for a single seed on each read
    // SeedId(tile, cluster, seed, reverse)
    // ReferencePosition(contigId, position)
    using isaac::alignment::SeedId;
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, s0, false), ReferencePosition(4, 20))));
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, s1, true ), ReferencePosition(4, 6))));
    using isaac::alignment::FragmentBuilder;
    using isaac::alignment::Cigar;
    FragmentBuilder fragmentBuilder(flowcells, 123, seedMetadataList.size()/2, 8, false,
                                    ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE,
                                    ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    // build the fragments
    fragmentBuilder.build(contigList, readMetadataList, seedMetadataList, testAdapters, matchList.begin(), matchList.end(), cluster4l, true);
    // Fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((long)0, fragmentBuilder.getFragments()[0][0].position);
    CPPUNIT_ASSERT_EQUAL(56U, fragmentBuilder.getFragments()[0][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragmentBuilder.getFragments()[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragmentBuilder.getFragments()[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragmentBuilder.getFragments()[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((44<<4)|Cigar::SOFT_CLIP), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((56<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].mismatchCount);
//    CPPUNIT_ASSERT_EQUAL((unsigned)44, fragmentBuilder.getFragments()[0][0].mismatchCount);
    // Fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((long)2, fragmentBuilder.getFragments()[1][0].position);
    CPPUNIT_ASSERT_EQUAL(58U, fragmentBuilder.getFragments()[1][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragmentBuilder.getFragments()[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragmentBuilder.getFragments()[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragmentBuilder.getFragments()[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragmentBuilder.getFragments()[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((58<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[2]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((42<<4)|Cigar::SOFT_CLIP), fragmentBuilder.getCigarBuffer()[3]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[1][0].mismatchCount);
//    CPPUNIT_ASSERT_EQUAL((unsigned)42, fragmentBuilder.getFragments()[1][0].mismatchCount);
}

void TestFragmentBuilder::testTrailingSoftClips()
{
    const unsigned s0 = 0;
    const unsigned s1 = 3;
    // build fragments for a single seed on each read
    // SeedId(tile, cluster, seed, reverse)
    // ReferencePosition(contigId, position)
    using isaac::alignment::SeedId;
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, s0, false), ReferencePosition(4, 16))));
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, s1, true ), ReferencePosition(4, 10))));
    using isaac::alignment::FragmentBuilder;
    using isaac::alignment::Cigar;
    FragmentBuilder fragmentBuilder(flowcells, 123, seedMetadataList.size()/2, 8, false,
                                    ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE,
                                    ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    // build the fragments
    fragmentBuilder.build(contigList, readMetadataList, seedMetadataList, testAdapters, matchList.begin(), matchList.end(), cluster4t, true);
    // Fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((long)16, fragmentBuilder.getFragments()[0][0].position);
    CPPUNIT_ASSERT_EQUAL(44U, fragmentBuilder.getFragments()[0][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragmentBuilder.getFragments()[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragmentBuilder.getFragments()[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragmentBuilder.getFragments()[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((44<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((56<<4)|Cigar::SOFT_CLIP), fragmentBuilder.getCigarBuffer()[1]);
//    CPPUNIT_ASSERT_EQUAL((unsigned)56, fragmentBuilder.getFragments()[0][0].mismatchCount);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].mismatchCount);
    // Fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((long)0, fragmentBuilder.getFragments()[1][0].position);
    CPPUNIT_ASSERT_EQUAL(42U, fragmentBuilder.getFragments()[1][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragmentBuilder.getFragments()[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragmentBuilder.getFragments()[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragmentBuilder.getFragments()[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragmentBuilder.getFragments()[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((58<<4)|Cigar::SOFT_CLIP), fragmentBuilder.getCigarBuffer()[2]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((42<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[3]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[1][0].mismatchCount);
//    CPPUNIT_ASSERT_EQUAL((unsigned)58, fragmentBuilder.getFragments()[1][0].mismatchCount);
}

void TestFragmentBuilder::testLeadingAndTrailingSoftClips()
{
    const unsigned s0 = 1;
    const unsigned s1 = 4;
    // build fragments for a single seed on each read
    // SeedId(tile, cluster, seed, reverse)
    // ReferencePosition(contigId, position)
    using isaac::alignment::SeedId;
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, s0, false), ReferencePosition(4, 22))));
    matchList.push_back(Match(Match(SeedId(tile0, 0, clusterId0, s1, true ), ReferencePosition(4, 11))));
    using isaac::alignment::FragmentBuilder;
    using isaac::alignment::Cigar;
    FragmentBuilder fragmentBuilder(flowcells, 123, seedMetadataList.size()/2, 8, false,
                                    ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE,
                                    ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    // build the fragments
    fragmentBuilder.build(contigList, readMetadataList, seedMetadataList, testAdapters, matchList.begin(), matchList.end(), cluster4lt, true);
    // Fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((long)0, fragmentBuilder.getFragments()[0][0].position);
    CPPUNIT_ASSERT_EQUAL(60U, fragmentBuilder.getFragments()[0][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragmentBuilder.getFragments()[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragmentBuilder.getFragments()[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragmentBuilder.getFragments()[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((10<<4)|Cigar::SOFT_CLIP), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((60<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((30<<4)|Cigar::SOFT_CLIP), fragmentBuilder.getCigarBuffer()[2]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[0][0].mismatchCount);
//    CPPUNIT_ASSERT_EQUAL((unsigned)40, fragmentBuilder.getFragments()[0][0].mismatchCount);
    // Fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((long)0, fragmentBuilder.getFragments()[1][0].position);
    CPPUNIT_ASSERT_EQUAL(60U, fragmentBuilder.getFragments()[1][0].observedLength);
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragmentBuilder.getFragments()[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragmentBuilder.getFragments()[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragmentBuilder.getFragments()[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragmentBuilder.getFragments()[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), fragmentBuilder.getCigarBuffer()[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((25<<4)|Cigar::SOFT_CLIP), fragmentBuilder.getCigarBuffer()[3]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((60<<4)|Cigar::ALIGN), fragmentBuilder.getCigarBuffer()[4]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((15<<4)|Cigar::SOFT_CLIP), fragmentBuilder.getCigarBuffer()[5]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragmentBuilder.getFragments()[1][0].mismatchCount);
//    CPPUNIT_ASSERT_EQUAL((unsigned)40, fragmentBuilder.getFragments()[1][0].mismatchCount);
}

