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
 **
 ** \file testSimpleIndelAligner.cpp
 **
 ** Tests for medium-size gap aligner.
 **
 ** \author Roman Petrovski
 **/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testSimpleIndelAligner.hh"

#include "alignment/fragmentBuilder/SimpleIndelAligner.hh"
#include "alignment/Cluster.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestSimpleIndelAligner, registryName("SimpleIndelAligner"));

static isaac::flowcell::ReadMetadataList getReadMetadataList(const unsigned readLength)
{
    std::vector<isaac::flowcell::ReadMetadata> ret =
        boost::assign::list_of
            (isaac::flowcell::ReadMetadata(1, readLength, 0, 0))
            ;
    return ret;
}

static isaac::alignment::SeedMetadataList getSeedMetadataList(const unsigned readLength)
{
    std::vector<isaac::alignment::SeedMetadata> ret =
        boost::assign::list_of
        (isaac::alignment::SeedMetadata( 0, 32, 0, 0))
        (isaac::alignment::SeedMetadata(readLength - 32 - 1, 32, 0, 1))
        ;
    return ret;
}

TestSimpleIndelAligner::TestSimpleIndelAligner()
{

}

void TestSimpleIndelAligner::setUp()
{
}

void TestSimpleIndelAligner::tearDown()
{
}

inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

namespace testSimpleIndelAligner
{

static const std::string irrelevantQualities("CFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDB<F>FGGBFGGFGCGGGDGGDDFHHHFEGGBGDGGBGGBEGEGGBGEHDHHHGGGGGDGGGG?GGGGCFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDBCFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDB");

struct ReadInit : public std::pair<std::string, std::string>
{
    static std::string reverseString(std::string fwd)
    {
        std::reverse(fwd.begin(), fwd.end());
        return fwd;
    }

    typedef std::pair<std::string, std::string> BaseType;
    ReadInit(const std::string &read, const bool reverse) :
        BaseType(reverse ? reverseString(read) : read, irrelevantQualities)
    {
        ISAAC_ASSERT_MSG(irrelevantQualities.length() >= read.length(), "Size matters");
    }
};

} // namespace TestSimpleIndelAligner

namespace isaac
{
namespace alignment
{

inline void phredToBcl(char &qual)
{
    qual -= 33;
}

template<class InpuT> InpuT& operator >>(InpuT &input, isaac::alignment::Read &read);

template<> testSimpleIndelAligner::ReadInit& operator >><testSimpleIndelAligner::ReadInit >(
    testSimpleIndelAligner::ReadInit &input,
    isaac::alignment::Read &read)
{
    ISAAC_ASSERT_MSG(input.first.length() <= input.second.length(), "sequence and quality must be of equal lengths");

    read.forwardSequence_ = vectorFromString(input.first);
    read.forwardQuality_ = vectorFromString(input.second);
    read.forwardQuality_.resize(read.forwardSequence_.size());

    std::for_each(read.forwardQuality_.begin(), read.forwardQuality_.end(), &phredToBcl);

    read.reverseSequence_ = read.forwardSequence_;
    read.reverseQuality_ = read.forwardQuality_;
    std::reverse(read.reverseSequence_.begin(), read.reverseSequence_.end());
    std::reverse(read.reverseQuality_.begin(), read.reverseQuality_.end());

    return input;
}

}
}

static const isaac::reference::Contig makeContig(const std::string forward)
{
    isaac::reference::Contig ret(0, "vasja");
    ret.forward_ = vectorFromString(forward);
    return ret;
}

static const int MATCH_SCORE = 0;
static const int MISMATCH_SCORE = -1;
static const int GAP_OPEN_SCORE = -2;
static const int GAP_EXTEND_SCORE = -1;
static const int MIN_GAP_EXTEND_SCORE = -5;

void TestSimpleIndelAligner::align(
    const std::string &read,
    const std::string &reference,
    isaac::alignment::FragmentMetadataList &fragmentMetadataList)
{
    const long long pos = read.find_first_not_of(' ');
    const std::string readWithoutSpaces = read.substr(pos);
    align(read, reference, getSeedMetadataList(std::min<long>(readWithoutSpaces.length(), reference.length() - pos)), fragmentMetadataList);
}

class TestAligner : public isaac::alignment::fragmentBuilder::SimpleIndelAligner
{
public:
    TestAligner() : isaac::alignment::fragmentBuilder::SimpleIndelAligner(
        MATCH_SCORE, MISMATCH_SCORE, GAP_OPEN_SCORE, GAP_EXTEND_SCORE, MIN_GAP_EXTEND_SCORE, 20000){}


    unsigned updateFragmentCigar(
        const isaac::flowcell::ReadMetadataList &readMetadataList,
        const std::vector<char> &reference,
        isaac::alignment::FragmentMetadata &fragmentMetadata,
        long strandPosition,
        isaac::alignment::Cigar &cigarBuffer,
        const unsigned cigarOffset) const
    {
        return isaac::alignment::fragmentBuilder::SimpleIndelAligner::updateFragmentCigar(
            readMetadataList, reference, fragmentMetadata, strandPosition, cigarBuffer, cigarOffset);
    }
};

void TestSimpleIndelAligner::align(
    const std::string &read,
    const std::string &reference,
    const isaac::alignment::SeedMetadataList &seedMetadataList,
    isaac::alignment::FragmentMetadataList &fragmentMetadataList)
{

    const long long referenceOffset = reference.find_first_not_of(' ');
    const std::string referenceWithoutSpaces = reference.substr(referenceOffset);
    const long long pos = read.find_first_not_of(' ');
    const std::string readWithoutSpaces = read.substr(pos);

    //cluster is referenced from fragmentMetadataList by pointer. Don't destroy it.
    static isaac::alignment::Cluster cluster(1000);
    testSimpleIndelAligner::ReadInit init(readWithoutSpaces, false);
    init >> cluster.at(0);

    static isaac::flowcell::ReadMetadataList readMetadatList;
    readMetadatList = getReadMetadataList(cluster.at(0).getLength());

    fragmentMetadataList.resize(2);

    fragmentMetadataList[0].readIndex = 0;
    fragmentMetadataList[0].contigId = 0;
    fragmentMetadataList[0].position = pos - referenceOffset;
    fragmentMetadataList[0].firstSeedIndex = seedMetadataList.at(0).getIndex();

    fragmentMetadataList[1].readIndex = 0;
    fragmentMetadataList[1].contigId = 0;
    fragmentMetadataList[1].position = reference.length() - readWithoutSpaces.length() - referenceOffset;
    fragmentMetadataList[1].firstSeedIndex = seedMetadataList.at(1).getIndex();

    const TestAligner aligner;

    const std::vector<char> referenceV(referenceWithoutSpaces.begin(), referenceWithoutSpaces.end());

    BOOST_FOREACH(isaac::alignment::FragmentMetadata &fragmentMetadata, fragmentMetadataList)
    {
        fragmentMetadata.cluster = &cluster;
        fragmentMetadata.cigarBuffer = &cigarBuffer_;
        fragmentMetadata.cigarOffset = cigarBuffer_.size();
        fragmentMetadata.observedLength = cluster.at(0).getLength();

        if (0 > fragmentMetadata.position)
        {
            cigarBuffer_.addOperation(std::max<long>(-fragmentMetadata.position, fragmentMetadata.leftClipped()), isaac::alignment::Cigar::SOFT_CLIP);
            ++fragmentMetadata.cigarLength;
            fragmentMetadata.observedLength -= std::max<long>(-fragmentMetadata.position, fragmentMetadata.leftClipped());
            fragmentMetadata.position += std::max<long>(-fragmentMetadata.position, fragmentMetadata.leftClipped());
        }
        else if (fragmentMetadata.leftClipped())
        {
            cigarBuffer_.addOperation(fragmentMetadata.leftClipped(), isaac::alignment::Cigar::SOFT_CLIP);
            ++fragmentMetadata.cigarLength;
            fragmentMetadata.observedLength -= fragmentMetadata.leftClipped();
            fragmentMetadata.position += fragmentMetadata.leftClipped();
        }

        long rightClip = 0;
        if (fragmentMetadata.rightClipped() || (fragmentMetadata.position + fragmentMetadata.observedLength > long(referenceV.size())))
        {
            rightClip = std::max<long>(fragmentMetadata.rightClipped(),
                                       (fragmentMetadata.position + fragmentMetadata.observedLength - referenceV.size()));
            fragmentMetadata.observedLength -= rightClip;
        }
        cigarBuffer_.addOperation(fragmentMetadata.observedLength, isaac::alignment::Cigar::ALIGN);
        ++fragmentMetadata.cigarLength;
        if (rightClip)
        {
            cigarBuffer_.addOperation(rightClip, isaac::alignment::Cigar::SOFT_CLIP);
            ++fragmentMetadata.cigarLength;
        }
        aligner.updateFragmentCigar(readMetadatList, referenceV, fragmentMetadata, fragmentMetadata.position, cigarBuffer_, fragmentMetadata.cigarOffset);
    }

    if (fragmentMetadataList[1].getUnclippedPosition() < fragmentMetadataList[0].getUnclippedPosition())
    {
        using std::swap;
        swap(fragmentMetadataList[1], fragmentMetadataList[0]);
    }

    isaac::reference::ContigList contigList;
    contigList.push_back(makeContig(referenceWithoutSpaces));

    aligner.alignSimpleIndels(cigarBuffer_, contigList,
                              readMetadatList,
                              seedMetadataList, fragmentMetadataList);
}

void TestSimpleIndelAligner::testEverything()
{
    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList;

        align("ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTAT"
                                                                                                   "AAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAAAAAAAAAAAAAAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("71M14D68M"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[0].getEditDistance());
    }

    { // verify proper preservation of alignment-independent clipping
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].leftClipped() = 8;
        fragmentMetadataList[1].rightClipped() = 7;

        std::vector<isaac::alignment::SeedMetadata> seedMetadataList =
            boost::assign::list_of
            // place seeds outside of clipped flanks
            (isaac::alignment::SeedMetadata( 10, 32, 0, 0))
            (isaac::alignment::SeedMetadata( 97, 32, 0, 1))
            ;

        align("ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTAT"
                                                                                                   "AAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAAAAAAAAAAAAAAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              seedMetadataList, fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("8S63M14D61M7S"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[0].getEditDistance());
        CPPUNIT_ASSERT_EQUAL(8U, unsigned(fragmentMetadataList[0].leftClipped()));
        CPPUNIT_ASSERT_EQUAL(7U, unsigned(fragmentMetadataList[0].rightClipped()));
    }


    {   // verify that despite the first 32 bases are occupied by a matching seed, the earliest possible deletion
        // position is selected at offset 25
        isaac::alignment::FragmentMetadataList fragmentMetadataList;

        align("GGTGCAGACTAGTAACAGTTGGTGGGCCGGCA"
                                                                                 "CTGATCCATTAATATATATTGCCCAGGTGCCGTGGCTCACCTATAATCCCAGCACTTTGAGAGGCCAA",
              "GGTGCAGACTAGTAACAGTTGGTGGGCCGGCACTGATCCATTAATATATATTGCCCAGGTGCCGGCACTGATCCATTAATATATATTGCCCAGGTGCCGTGGCTCACCTATAATCCCAGCACTTTGAGAGGCCAA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("25M35D75M"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(35U, fragmentMetadataList[0].getEditDistance());
    }

    {   // verify that despite the first 32 bases are occupied by a matching seed, the earliest possible deletion
        // position is selected at offset 0 and causes alignment position change instead of cigar beginning with a deletion
        isaac::alignment::FragmentMetadataList fragmentMetadataList;

        align("GGTGCAGACTAGTAACAGTTGGTGGGCCGGCA"
                                                                              "CTGATCCATTAATATATATTGCCCAGGTGCCGTGGCTCACCTATAATCCCAGCACTTTGAGAGGCCAA",
              "GGTGCAGACTAGTAACAGTTGGTGGGCCGGCAGGTGCAGACTAGTAACAGTTGGTGGGCCGGCACTGATCCATTAATATATATTGCCCAGGTGCCGTGGCTCACCTATAATCCCAGCACTTTGAGAGGCCAA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("100M"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(32UL, fragmentMetadataList[0].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[0].getEditDistance());
    }

    { //verifying that it picks the earliest possible position for the insertion gap
        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        std::vector<isaac::alignment::SeedMetadata> seedMetadataList =
            boost::assign::list_of
            // place seeds outside of clipped flanks
            (isaac::alignment::SeedMetadata( 0, 32, 0, 0))
            (isaac::alignment::SeedMetadata( 87, 32, 0, 1))
            ;
        align("ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTAT"
                                                                       "AAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              seedMetadataList, fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("57M14I68M"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[0].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        std::vector<isaac::alignment::SeedMetadata> seedMetadataList =
            boost::assign::list_of
            // place seeds outside of clipped flanks
            (isaac::alignment::SeedMetadata( 0, 32, 0, 0))
            (isaac::alignment::SeedMetadata( 87, 32, 0, 1))
            ;

        align("ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTAT"
                                                                       "TAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACTAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              seedMetadataList, fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("57M14I68M"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[0].getEditDistance());
    }

    { // verify proper preservation of alignment-independent clipping for insertions
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].leftClipped() = 8;
        fragmentMetadataList[1].rightClipped() = 7;

        std::vector<isaac::alignment::SeedMetadata> seedMetadataList =
            boost::assign::list_of
            // place seeds outside of clipped flanks
            (isaac::alignment::SeedMetadata( 10, 32, 0, 0))
            (isaac::alignment::SeedMetadata( 87, 32, 0, 1))
            ;

        align("ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTAT"
                                                                       "TAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACTAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              seedMetadataList, fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("8S49M14I61M7S"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[0].getEditDistance());
        CPPUNIT_ASSERT_EQUAL(8U, unsigned(fragmentMetadataList[0].leftClipped()));
        CPPUNIT_ASSERT_EQUAL(7U, unsigned(fragmentMetadataList[0].rightClipped()));
    }

    { //verifying that it picks the earliest possible position for the insertion but not within the anchoring seed
        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align("TTCTGTCTCTTATACAACAAGTGGATGTGTAA"
                                "AAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "TTCTGTCTCTTATACAACAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("32M14I54M"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(8U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(22U, fragmentMetadataList[0].getEditDistance());
    }

    {

        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align("              AAAAAAAAAAAAAATTCTGTCTCTTATACAAC"
                                              "AAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "TTTTTTTTTTTTTTAAAAAAAAAAAAAATTCTGTCTCTTATACAACCCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("32M14I54M"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(14UL, fragmentMetadataList[0].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[0].getEditDistance());
    }

    { // ensure that the insertion is not placed before head seed or on tail seed
        std::vector<isaac::alignment::SeedMetadata> seedMetadataList =
            boost::assign::list_of
            (isaac::alignment::SeedMetadata( 64, 32, 0, 0))
            (isaac::alignment::SeedMetadata( 32, 32, 0, 1))
            ;

        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align(                             "TGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCT"
                                                        "CTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              seedMetadataList,
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("128M"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0UL, fragmentMetadataList[0].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(39U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(39U, fragmentMetadataList[0].getEditDistance());
    }

    { // ensure that the insertion is not placed before head seed
        std::vector<isaac::alignment::SeedMetadata> seedMetadataList =
            boost::assign::list_of
            (isaac::alignment::SeedMetadata( 64, 32, 0, 0))
            (isaac::alignment::SeedMetadata( 0, 32, 0, 1))
            ;

        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align(                             "TGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCT"
                                                        "CTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              seedMetadataList,
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("32M29I67M"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(29UL, fragmentMetadataList[0].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(21U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(50U, fragmentMetadataList[0].getEditDistance());
    }
    
    { // ensure that reference start soft clipping does not break
        std::vector<isaac::alignment::SeedMetadata> seedMetadataList =
            boost::assign::list_of
            (isaac::alignment::SeedMetadata( 32, 32, 0, 0))
            (isaac::alignment::SeedMetadata( 64, 32, 0, 1))
            ;

        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTC"
                                                                                 "TCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              "                     GGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCAAATCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              seedMetadataList,
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("21S43M3D93M"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0UL, fragmentMetadataList[0].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(3U, fragmentMetadataList[0].getEditDistance());
    }

    { // ensure that reference start soft clipping does not break
        std::vector<isaac::alignment::SeedMetadata> seedMetadataList =
            boost::assign::list_of
            (isaac::alignment::SeedMetadata( 32, 32, 0, 0))
            (isaac::alignment::SeedMetadata( 67, 32, 0, 1))
            ;

        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTC"
                                                                           "TCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              "                     GGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              seedMetadataList,
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("21S43M3I90M"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0UL, fragmentMetadataList[0].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(3U, fragmentMetadataList[0].getEditDistance());
    }

    { // ensure that reference start soft clipping does not break
        std::vector<isaac::alignment::SeedMetadata> seedMetadataList =
            boost::assign::list_of
            (isaac::alignment::SeedMetadata( 28, 32, 0, 0))
            (isaac::alignment::SeedMetadata( 64, 32, 0, 1))
            ;

        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGACTC"
                                                                           "TCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              "                 TCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGATCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              seedMetadataList,
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::string("17S44M3I93M"), fragmentMetadataList[0].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0UL, fragmentMetadataList[0].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[0].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(3U, fragmentMetadataList[0].getEditDistance());
    }


    { // avoid gaps with too many mismatches around them

        { // make sure gap gets accepted when there are less than 8 mismatches on each side
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACTTTCCCCTCTTGGGCCTTCTGGGAACGACCCCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCT"
                   /*                 x           xxxxxxxx x                                */     "GCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::string("73M8D76M"), fragmentMetadataList[0].getCigarString());
            CPPUNIT_ASSERT_EQUAL(0UL, fragmentMetadataList[0].getFStrandReferencePosition().getPosition());
            CPPUNIT_ASSERT_EQUAL(11U, fragmentMetadataList[0].getMismatchCount());
            CPPUNIT_ASSERT_EQUAL(19U, fragmentMetadataList[0].getEditDistance());
        }

        { // too many mismatches in the right flank of the gap
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACTTTCCCCTCTTGGGCCTTCTGGGAACGACCCCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCT"     /*xxxxxxx                xx*/
                   /*                 x           xxxxxxxx x                                */     "ATGTGATGCCTCTCTGCGCCTGCGTCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::string("149M"), fragmentMetadataList[0].getCigarString());
        }

        { // too many mismatches in the left flank of the gap
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACTTTCCCCTCTTGGGCCTTCTGGGAACGACCCCCTCCGCTGGGGCGGAGGTCCTCACCGCGACT"
                   /*                 x           xxxxxxxx x   x   x  x  x x  x    x x   x  */     "TGTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::string("149M"), fragmentMetadataList[0].getCigarString());
        }

        { // ok mismatches in the right flank of the gap
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACTTTCCCCTCTTGGGCCTTCTGGGAACGACCCCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCT"     /*xxx xxx                xx*/
                   /*                 x           xxxxxxxx x                                */     "ATGGGATGCCTCTCTGCGCCTGCGTCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::string("73M8D76M"), fragmentMetadataList[0].getCigarString());
        }

        { // ok mismatches in the left flank of the gap
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACTTTCCCCTCTTGGGCCTTCTGGGAACGACCCCCTCCGCTGGGGCGGAGGTCCTCTCCGCGACT"
                   /*                 x           xxxxxxxx x   x   x  x  x x  x      x   x  */     "TGTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::string("73M8D76M"), fragmentMetadataList[0].getCigarString());
        }
    }

    { // avoid gaps with too many mismatches around them

        { // make sure gap gets accepted when there are 8 mismatches on each side
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCCCCGTCGGGCCCTTTCCCCCTGCCC"/*x x x x x x x x*/
                  /*                                                   x x  x  x x  x x x */"AAACGGGGGTCTCGCCGTGTGTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::string("77M3I80M"), fragmentMetadataList[0].getCigarString());
            CPPUNIT_ASSERT_EQUAL(0UL, fragmentMetadataList[0].getFStrandReferencePosition().getPosition());
            CPPUNIT_ASSERT_EQUAL(16U, fragmentMetadataList[0].getMismatchCount());
            CPPUNIT_ASSERT_EQUAL(19U, fragmentMetadataList[0].getEditDistance());
        }

        { // too many mismatches in the right flank of the gap
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCCCCGTCGGGCCCTTTCCCCCTGCCC"/*x x x x x x x x x*/
                  /*                                                   x x  x  x x  x x x */"AAACGGGGGTCTCGCCGTGTGTCCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::string("3S157M"), fragmentMetadataList[0].getCigarString());
        }

        { // too many mismatches in the left flank of the gap
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCCCCGTCGGGCCTTTTCCCCCTGCCC"/*x x x x x x x x*/
                  /*                                                   x x  x  x xx x x x */"AAACGGGGGTCTCGCCGTGTGTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::string("3S157M"), fragmentMetadataList[0].getCigarString());
        }

    }

}

