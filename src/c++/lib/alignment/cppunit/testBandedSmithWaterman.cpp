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

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/assign.hpp>
#include <boost/assign/std/vector.hpp> 

using namespace std;

#include "RegistryName.hh"
#include "testBandedSmithWaterman.hh"
#include "BuilderInit.hh"
#include "alignment/Cigar.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestBandedSmithWaterman, registryName("BandedSmithWaterman"));

std::string getGenome(unsigned size = 1000)
{
    static const std::string bases = "ACGT";
    std::string genome;
    genome.reserve(size);
    while (size--)
    {
        genome.push_back(bases[rand() % bases.size()]);
    }
    return genome;
}

TestBandedSmithWaterman::TestBandedSmithWaterman()
    : bsw(2, -1, 15, 3, 300)
//    : bsw(2, -1, 15, 3, 300)
    , genome(getGenome())
{
}

void TestBandedSmithWaterman::setUp()
{
}

void TestBandedSmithWaterman::tearDown()
{
}


void TestBandedSmithWaterman::testCustom()
{
    //const std::string query = "GGCAATAACCTAAGATAAAAAAAAACATTGTAAATAAAAGTAAATGCTGTAATATTTAAGGTAAGTAAGCCAACAGCTGAGGAAAAGGGAAATCCATTAT";
    //const std::string database = "CTTTAATCAGGCAATAACCTAAGATAAAAAAAACATTGTAAACAAAAGTAAATGCTGTAATATTTAAGGTAAGTAAGCCAACAGCTGAGGAAAAGGGAAATCCATTATCCTGGGT";
    const std::vector<char> query = vectorFromString("CTAAGACCCCACACTCTGGGACACCAAGGTGGGAGGATCGCTGGAGCTCAGGAGTTTGAGACCAGCCTGGACAACATGGTGTGACCCTGTCTACAGAAAA");
    const std::vector<char> database = vectorFromString("AATGCCTCTGGCCTGGGCGTGGGAGTTCATGCTTGTAATCGCATATCGCTAGAGCCCAGGAGTTTGAGACCAGCCTGGACAACATGGTGAAAACCCTCGTTGCTACTAAAAATAC");
    //ii = 27
    isaac::alignment::Cigar cigar;
    cigar.clear();
    bsw.align(query, database.begin(), database.end(), cigar);
    CPPUNIT_ASSERT_EQUAL(5UL, cigar.size());
    CPPUNIT_ASSERT_EQUAL((9 << 4) + 2U, cigar[0]);
    CPPUNIT_ASSERT_EQUAL((16 << 4) + 0U, cigar[1]);
    CPPUNIT_ASSERT_EQUAL((1 << 4) + 1U, cigar[2]);
    CPPUNIT_ASSERT_EQUAL((83 << 4) + 0U, cigar[3]);
    CPPUNIT_ASSERT_EQUAL((7 << 4) + 2U, cigar[4]);
}

void TestBandedSmithWaterman::testUngapped()
{
    const std::vector<char> database = subv(genome, 100, 115);
    std::vector<char> query = subv(database, 0, 100);
    CPPUNIT_ASSERT_EQUAL(115UL, database.size());
    isaac::alignment::Cigar cigar;
    bsw.align(query, database.begin(), database.end(), cigar);
    CPPUNIT_ASSERT_EQUAL(1UL, cigar.size());
    CPPUNIT_ASSERT_EQUAL(1600U, cigar[0]);
    CPPUNIT_ASSERT_EQUAL(242U, cigar[1]);
    cigar.clear();
    for (unsigned i = 1; i < 15; ++i)
    {
        query = subv(database, i, 100);
        bsw.align(query, database.begin(), database.end(), cigar);
        CPPUNIT_ASSERT_EQUAL(1UL, cigar.size());
        CPPUNIT_ASSERT_EQUAL(1600U, cigar[0]);
        cigar.clear();
    }
    query = subv(database, 15, 100);
    bsw.align(query, database.begin(), database.end(), cigar);
    CPPUNIT_ASSERT_EQUAL(1UL, cigar.size());
    CPPUNIT_ASSERT_EQUAL(1600U, cigar[0]);
    cigar.clear();
}

void TestBandedSmithWaterman::testSingleDeletion()
{
    // Note: Deletions in homopolymer won't show at the right place
    // Note: the gap penalties make it impractical to have inserts greater than 9
    const unsigned left = 40;
    const unsigned right = 40;
    const std::string deletion = "AGAGCAGCGAGCGACAGCAGCAGCAAA"; // no Ts to ensure constant location
    for (unsigned deletionLength = 1; 13 >= deletionLength; ++deletionLength)
    {
        const unsigned dl = 7 - (deletionLength / 2);
        const std::vector<char> dlS = subv(genome, 100, dl);
        const std::vector<char> leftS = subv(genome, 100 + dl, left - 1) + "T";
        const std::vector<char> rightS = subv(genome, 100 + dl + left, right);
        const std::vector<char> deletionS1 = subv(deletion, 0, deletionLength);
        const unsigned drD = 15 - dl - deletionS1.size();
        const std::vector<char> drDS = subv(genome, 100 + dl + left + right, drD);
        const std::vector<char> databaseD = dlS + leftS + deletionS1 + rightS + drDS;
        const std::vector<char> queryD = leftS + rightS;
        isaac::alignment::Cigar cigar;
        bsw.align(queryD, databaseD.begin(), databaseD.end(), cigar);
        CPPUNIT_ASSERT_EQUAL(3UL, cigar.size());
        CPPUNIT_ASSERT_EQUAL(( left << 4) | 0U, cigar[0]);
        CPPUNIT_ASSERT_EQUAL((unsigned)( deletionS1.size() << 4) | 2U, cigar[1]);
        CPPUNIT_ASSERT_EQUAL((right << 4) | 0U, cigar[2]);
        cigar.clear();
    }
}

void TestBandedSmithWaterman::testSingleInsertion()
{
    // Note: Insertions require appropriate database length at the beginning
    const std::vector<char> database = subv(genome, 100, 220);
    unsigned queryLength = database.size() - 15;
    // Note: the gap penalties make it impractical to have inserts greater than 9
    for (unsigned insertLength = 1; 9 >= insertLength; ++insertLength)
    {
        unsigned left = 100;
        unsigned right = queryLength - left - insertLength;
        unsigned dl = 9;
        std::vector<char> query = subv(database, dl, left) + std::string(insertLength, 'T') + subv(database, left + dl, right);
        //CPPUNIT_ASSERT_EQUAL(115UL, database.size());
        isaac::alignment::Cigar cigar;
        bsw.align(query, database.begin(), database.end(), cigar);
        CPPUNIT_ASSERT_EQUAL(3UL, cigar.size());
        CPPUNIT_ASSERT_EQUAL(( left << 4) | 0U, cigar[0]);
        CPPUNIT_ASSERT_EQUAL(( insertLength << 4) | 1U, cigar[1]);
        CPPUNIT_ASSERT_EQUAL((right << 4) | 0U, cigar[2]);
        cigar.clear();
    }
}

void TestBandedSmithWaterman::testMultipleIndels()
{
    // Note: the gap penalties make it impractical to have inserts greater than 9
    const unsigned left = 20;
    const unsigned center = 20;
    const unsigned right = 20;
    const unsigned dl = 6;
    const std::string dlS = genome.substr(100, dl);
    const std::string leftS = genome.substr(100 + dl, left - 1) + "T";
    const std::string insertS1 = "A";
    const std::string insertS2 = "CG";
    const std::string centerS = genome.substr(100 + dl + left, center - 1) + "T";
    const std::string deletionS1 = "AAG";
    const std::string deletionS2 = "ACAG";
    const std::string rightS = genome.substr(100 + dl + left + center, right);
    // 1 Insertion and 1 Deletion
    const unsigned drID = 15 - dl + insertS1.length() - deletionS2.length();
    const std::string drIDS = genome.substr(100 + dl + left + center + right, drID);
    const std::vector<char> databaseID = vectorFromString(dlS + leftS + centerS + deletionS2 + rightS + drIDS);
    const std::vector<char> queryID = vectorFromString(leftS + insertS1 + centerS + rightS);
    //CPPUNIT_ASSERT_EQUAL(115UL, database.size());
    isaac::alignment::Cigar cigar;
    bsw.align(queryID, databaseID.begin(), databaseID.end(), cigar);
    CPPUNIT_ASSERT_EQUAL(5UL, cigar.size());
    CPPUNIT_ASSERT_EQUAL(( left << 4) | 0U, cigar[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)( insertS1.length() << 4) | 1U, cigar[1]);
    CPPUNIT_ASSERT_EQUAL((center << 4) | 0U, cigar[2]);
    CPPUNIT_ASSERT_EQUAL((unsigned)(deletionS2.length() << 4) | 2U, cigar[3]);
    CPPUNIT_ASSERT_EQUAL((right << 4) | 0U, cigar[4]);
    cigar.clear();
    // 2 Insertions
    const unsigned drI2 = 15 - dl + insertS1.length() +insertS2.length();
    const std::string drI2S = genome.substr(100 + dl + left + center + right, drI2);
    const std::vector<char> databaseI2 = vectorFromString(dlS + leftS + centerS + rightS + drI2S);
    const std::vector<char> queryI2 = vectorFromString(leftS + insertS1 + centerS + insertS2 + rightS);
    bsw.align(queryI2, databaseI2.begin(), databaseI2.end(), cigar);
    CPPUNIT_ASSERT_EQUAL(5UL, cigar.size());
    CPPUNIT_ASSERT_EQUAL(( left << 4) | 0U, cigar[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)(insertS1.length() << 4) | 1U, cigar[1]);
    CPPUNIT_ASSERT_EQUAL((center << 4) | 0U, cigar[2]);
    CPPUNIT_ASSERT_EQUAL((unsigned)(insertS2.length() << 4) | 1U, cigar[3]);
    CPPUNIT_ASSERT_EQUAL((right << 4) | 0U, cigar[4]);
    cigar.clear();
    // 2 deletions
    const unsigned drD2 = 15 - dl - deletionS1.length() - deletionS2.length();
    const std::string drD2S = genome.substr(100 + dl + left + center + right, drD2);
    const std::vector<char> databaseD2 = vectorFromString(dlS + leftS + deletionS1 + centerS + deletionS2 + rightS + drD2S);
    const std::vector<char> queryD2 = vectorFromString(leftS + centerS + rightS);
    bsw.align(queryD2, databaseD2.begin(), databaseD2.end(), cigar);
    CPPUNIT_ASSERT_EQUAL(5UL, cigar.size());
    CPPUNIT_ASSERT_EQUAL(( left << 4) | 0U, cigar[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)(deletionS1.length() << 4) | 2U, cigar[1]);
    CPPUNIT_ASSERT_EQUAL((center << 4) | 0U, cigar[2]);
    CPPUNIT_ASSERT_EQUAL((unsigned)(deletionS2.length() << 4) | 2U, cigar[3]);
    CPPUNIT_ASSERT_EQUAL((right << 4) | 0U, cigar[4]);
    cigar.clear();
}
