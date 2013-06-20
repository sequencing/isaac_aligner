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
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testSequencingAdapter.hh"

#include "alignment/Cluster.hh"
#include "alignment/matchSelector/FragmentSequencingAdapterClipper.hh"
#include "flowcell/SequencingAdapterMetadata.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestSequencingAdapter, registryName("SequencingAdapter"));

inline isaac::flowcell::ReadMetadataList getReadMetadataList(const unsigned l0 = 100, const unsigned l1 = 100)
{
    std::vector<isaac::flowcell::ReadMetadata> ret =
        boost::assign::list_of
            (isaac::flowcell::ReadMetadata(1, l0, 0, 0))
            (isaac::flowcell::ReadMetadata(l0 + 1, l0 + l1, 1, l0))
            ;
    return ret;
}

inline isaac::alignment::SeedMetadataList getSeedMetadataList()
{
    std::vector<isaac::alignment::SeedMetadata> ret =
        boost::assign::list_of
        (isaac::alignment::SeedMetadata( 0, 32, 0, 0))
        (isaac::alignment::SeedMetadata(32, 32, 0, 1))
        (isaac::alignment::SeedMetadata(64, 32, 0, 2))
        (isaac::alignment::SeedMetadata( 0, 32, 1, 3))
        (isaac::alignment::SeedMetadata(32, 32, 1, 4))
        (isaac::alignment::SeedMetadata(64, 32, 1, 5))
        ;
    return ret;
}

static const isaac::flowcell::SequencingAdapterMetadata testAdapterMetadataLeft(
    "CTGTCTCTTATACACATCT",
    false,
    strlen("CTGTCTCTTATACACATCT"));

static const isaac::flowcell::SequencingAdapterMetadata testAdapterMetadataRight(
    "AGATGTGTATAAGAGACAG",
    true,
    strlen("AGATGTGTATAAGAGACAG"));

static const isaac::alignment::matchSelector::SequencingAdapterList matePairAdapters =
    boost::assign::list_of(isaac::alignment::matchSelector::SequencingAdapter(testAdapterMetadataLeft))
                          (isaac::alignment::matchSelector::SequencingAdapter(testAdapterMetadataRight));

static const isaac::flowcell::SequencingAdapterMetadata standardAdapterMetadataLeft(
    "CTGTCTCTTATACACATCT",
    false,
    0);

static const isaac::flowcell::SequencingAdapterMetadata standardAdapterMetadataRight(
    "AGATGTGTATAAGAGACAG",
    true,
    0);

static const isaac::alignment::matchSelector::SequencingAdapterList standardAdapters =
    boost::assign::list_of(isaac::alignment::matchSelector::SequencingAdapter(standardAdapterMetadataLeft))
                          (isaac::alignment::matchSelector::SequencingAdapter(standardAdapterMetadataRight));

static const int ELAND_MATCH_SCORE = 2;
static const int ELAND_MISMATCH_SCORE = -1;
static const int ELAND_GAP_OPEN_SCORE = -15;
static const int ELAND_GAP_EXTEND_SCORE = -3;
static const int ELAND_MIN_GAP_EXTEND_SCORE = 25;

TestSequencingAdapter::TestSequencingAdapter() :
    readMetadataList(getReadMetadataList()),
    seedMetadataList(getSeedMetadataList()),
    flowcells(1, isaac::flowcell::Layout("", isaac::flowcell::Layout::Fastq, std::vector<unsigned>(),
                                         readMetadataList, seedMetadataList, "blah")),
    ungappedAligner_(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE),
    irrelevantQualities("CFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDB<F>FGGBFGGFGCGGGDGGDDFHHHFEGGBGDGGBGGBEGEGGBGEHDHHHGGGGGDGGGG?GGGGDBEDDEGEHHFHEGBHHDDDB<F>FGGBFGGFGCGGGDGGDDFHHHFEGGBGDGDBEDDEGEHHFHEGBHHDDDB<F>FGGBFGGFGCGGGDGGDDFHHHFEGGBGDG")
{

}

void TestSequencingAdapter::setUp()
{
}

void TestSequencingAdapter::tearDown()
{
}

inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

namespace isaac
{
namespace alignment
{

inline void phredToBcl(char &qual)
{
    qual -= 33;
}

template<class InpuT> InpuT& operator >>(InpuT &input, isaac::alignment::Read &read);

template<> std::pair<std::string, std::string>& operator >><std::pair<std::string, std::string> >(
    std::pair<std::string, std::string> &input,
    isaac::alignment::Read &read)
{
    ISAAC_ASSERT_MSG(input.first.length() == input.second.length(), "sequence and quality must be of equal lengths");

    read.forwardSequence_ = vectorFromString(input.first);
    read.forwardQuality_ = vectorFromString(input.second);

    std::for_each(read.forwardQuality_.begin(), read.forwardQuality_.end(), &phredToBcl);

    read.reverseSequence_ = read.forwardSequence_;
    read.reverseQuality_ = read.forwardQuality_;
    std::reverse(read.reverseSequence_.begin(), read.reverseSequence_.end());
    std::reverse(read.reverseQuality_.begin(), read.reverseQuality_.end());

    return input;
}

}
}

void TestSequencingAdapter::testEverything()
{
    testMp51M49S();
    testMp51S49M();
    testMp94M6S();
    testMp33S67M();
    testMp40M60S();
    testMp47S53M();
    testMp30S70M();
    testMp11S89M();
    testMp16S84M();
    testStd38M62S();
    testStd76S24MReverse();
    testStd36M114S();
    testStdBeforeSequence();
    testStdReverseAfterSequence();
    testStdReverseSequenceTooGood();

}

static const isaac::reference::Contig makeContig(const std::string forward)
{
    isaac::reference::Contig ret(0, "vasja");
    ret.forward_ = vectorFromString(forward);
    return ret;
}

static std::string reverse(std::string fwd)
{
    std::reverse(fwd.begin(), fwd.end());
    return fwd;
}

void TestSequencingAdapter::align(
    const std::string &read,
    const std::string &reference,
    const isaac::alignment::matchSelector::SequencingAdapterList &adapters,
    isaac::alignment::FragmentMetadata &fragmentMetadata)
{
    isaac::alignment::Cluster cluster(isaac::flowcell::getMaxReadLength(flowcells));
    std::pair<std::string, std::string> init(fragmentMetadata.reverse ? reverse(read) : read,
                                             irrelevantQualities.substr(0, read.length()));
    init >> cluster.at(0);

    fragmentMetadata.contigId = 0;
    fragmentMetadata.position = 0;
    fragmentMetadata.cluster = &cluster;
    fragmentMetadata.cigarBuffer = &cigarBuffer_;

    const isaac::reference::Contig referenceContig = makeContig(reference);

    isaac::alignment::matchSelector::FragmentSequencingAdapterClipper adapterClipper(adapters);
    adapterClipper.checkInitStrand(fragmentMetadata, referenceContig);
    ungappedAligner_.alignUngapped(fragmentMetadata, cigarBuffer_, readMetadataList, adapterClipper, referenceContig);
}

void TestSequencingAdapter::testMp51M49S()
{

    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
//                                                            CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG
//                                                            ||||||||||||||||||||||||||||||||||||||
    align("CGATTGTCTTTGCTGCCAATTTTAGCGTTGGCGTTAACGTCATGCTTAAGCCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGCTGCTACGCCA",
//         ||||||||||||||||||||||||||||||||||||||||||||||||||||xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|xxxx|
          "CGATTGTCTTTGCTGCCAATTTTAGCGTTGGCGTTAACGTCATGCTTAAGCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
          matePairAdapters,
          fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("51M49S"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(51U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getStrandReferencePosition());
}

/**
sequencing direction ->
alignment direction ->
xxxxxxxxxxxxxx|xxxxxxxxxxxxxxxxxx|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
GCACAGCATTACACCTGTCTCTTATACACATCTCTGGAATATGATACACCGCCGAGAAATCATCACCTTAACCTCTGATAATCGTCATATACCGGACAAG
              CTGTCTCTTATACACATCT-------------------
**/
void TestSequencingAdapter::testMp33S67M()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    //                   CTGTCTCTTATACACATCT-------------------
    //                   |||||||||||||||||||
    align("GCACAGCATTACACCTGTCTCTTATACACATCTCTGGAATATGATACACCGCCGAGAAATCATCACCTTAACCTCTGATAATCGTCATATACCGGACAAG",
         //xxxxxxxxxxxxxx|xxxxxxxxxxxxxxxxxx|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          "AAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAACTGGAATATGATACACCGCCGAGAAATCATCACCTTAACCTCTGATAATCGTCATATACCGGACAAG",
          matePairAdapters,
          fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("33S67M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(67U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 33U), fragmentMetadata.getStrandReferencePosition());

}

void TestSequencingAdapter::testMp51S49M()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    //                  CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG                        CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG
    //                  ||||||||||||||||||||||||||||||||||||||                                           ||||||
    align("CTCAGCCGTGAAGCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGCGCCCTACACCATAACCAGCCTGTAAAGAATAAGCGCCCATAAAGATGT",
         //xxxxxxxxxxxxx||xxxxxxxxxxxxxxxxx|x||||x|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          "AAAAAAAAAAAAACTAAAAAAAAAAAAAAAAAATATGTCTATAAGAGACAGCGCCCTACACCATAACCAGCCTGTAAAGAATAAGCGCCCATAAAGATGT",
          matePairAdapters,
          fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("51S49M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(49U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 51U), fragmentMetadata.getStrandReferencePosition());
}

void TestSequencingAdapter::testMp94M6S()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    //                  CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG                        CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG
    //                  |||||||||| ||||||||||||||||| |||||||||                                           ||||||
    align("CTCAGCCGTGAAGCTGTCTCTTAAACACATCTAGATGTGTAAAAGAGACAGCGCCCTACACCATAACCAGCCTGTAAAGAATAAGCGCCCATAAAGATGT",
         //xxxxxxxxxxxxx||xxxxxxxxxxxxxxxxx|x||||x|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          "AAAAAAAAAAAAACTAAAAAAAAAAAAAAAAAAGATGTGTATAAGAGACAGCGCCCTACACCATAACCAGCCTGTAAAGAATAAGCGCCCATAAACTATT",
          matePairAdapters,
          fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("94M6S"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(23U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(23U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(94U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getStrandReferencePosition());
}

void TestSequencingAdapter::testMp40M60S()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    //                                             CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG
    //                                             ||||||||||||||||||||||||||||||||||||||----------------------
    align("TGGGCCAGCTTCATGACATAACGCGGTTGTTGAGATAAAGCTGTCTCTTATACACATCTCTGACCAACCCAACGCCAGTCTTCGCCCCCTCCAGTTAACT",
         //|||||||||||||||||||||||||||||||||||||||||   |   |       |     |                      ||  | |      ||
          "TGGGCCAGCTTCATGACATAACGCGGTTGTTGAGATAAAGCGTCCAGCTTCGGCATTAATAAAGTTTGTGCGGCGTTATAAAAAACCGGTTCGAGATTCT",
          matePairAdapters,
          fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("40M60S"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(40U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getStrandReferencePosition());
}

void TestSequencingAdapter::testMp47S53M()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    //              CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG
    //     ---------||||||||||||||||||||||||||||||||||||||
    align("AGAGCTGGCCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTTCTTCACCCCCGCACCATTACCCCCATCGCCCAGTTCCAGATCCCTTGCCT",
         //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|||||||||||||||||||||||||||||||||||||||||||||||||||||
          "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTTCTTCACCCCCGCACCATTACCCCCATCGCCCAGTTCCAGATCCCTTGCCT",
          matePairAdapters,
          fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("47S53M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(53U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 47U), fragmentMetadata.getStrandReferencePosition());
}

void TestSequencingAdapter::testMp30S70M()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
// CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG
//         ||||||||||||||||||||||||||||||
    align("TATACACATCTAGATGTGTATAAGAGACAGGTGTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTCTCTTCTCTGGAATATGATACACCGCC",
         //x|x|x|x|xxx|x|xxxxx|x||x|x|x|x||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTCTCTTCTCTGGAATATGATACACCGCC",
          matePairAdapters,
          fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("30S70M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(70U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 30U), fragmentMetadata.getStrandReferencePosition());
}

void TestSequencingAdapter::testMp11S89M()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
// CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG
//         |||||||||||||||||||| |||||||||
    align("TATACACATCTAGATGTGTAAAAGAGACAGGTGTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTCTCTTCTCTGGAATATGATACACCGCC",
         // | | | |   | |     |||| | | | ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTCTCTTCTCTGGAATATGATACACCGCC",
          matePairAdapters,
          fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("11S89M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(10U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(10U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(89U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 11U), fragmentMetadata.getStrandReferencePosition());
}

void TestSequencingAdapter::testMp16S84M()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
//  CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG
//                        ||||||||||||||||
    align(               "TGTGTATAAGAGACAGGTGTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTCTCTTCTCTGGAATATGATACACCGCCTATACACATCTAGA",
                        //     | || | | | |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| | | | |   | |
                         "AAAAAAAAAAAAAAAAGTGTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTCTCTTCTCTGGAATATGATACACCGCCAAAAAAAAAAAAAA",
                         matePairAdapters,
                         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("16S84M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(8U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(8U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(84U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 16U), fragmentMetadata.getStrandReferencePosition());
}

void TestSequencingAdapter::testStd38M62S()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
                                             //  CTGTCTCTTATACACATCT*
                                             //                    *AGATGTGTATAAGAGACAG
                                             //  |||||||||||||||||||-------------------------------------------
    align("TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         //|||||||||||||||||||||||||||||||||||||||| |||  | |  |     | |     | || | | | |||||||||| | | | |   | |
          "TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTCTCTTCTCTGGAATATGATAAAAAAAAAAAAAAAAAGTGCACCGCCAAAAAAAAAAAAAA",
         standardAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("38M62S"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(38U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getStrandReferencePosition());
}

void TestSequencingAdapter::testStd76S24MReverse()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = true;
       //                                        CTGTCTCTTATACACATCT*
       //                                                          *AGATGTGTATAAGAGACAG
       //  ---------------------------------------------------------|||||||||||||||||||
    align("TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         //|| |||| ||| |||| ||| ||| | ||| |||| || | |||  | |  |     | |     | || | | | |||||||||| | | | |   | |
          "TGATTAATGTACCGGTTAAACCGTTTCACCTCAATTTTTTCTCTTCTCTGGAATATGATAAAAAAAAAAAAAAAAAGTGCACCGCCAAAAAAAAAAAAAA",
         standardAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("76S24M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(8U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(8U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(24U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 76U), fragmentMetadata.getFStrandReferencePosition());
}

void TestSequencingAdapter::testStd36M114S()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
                                             //CTGTCTCTTATACACATCT*
       //                                      |||||||||||||||||||-----------------------------------------------------------------------------------------------
    align("AGATAAGTCCATGAAGTCACCAGCACCGTCCATGTTCTGTCTCTTATACACATCTCCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAACACAACCATCGAGTCCACATCAGATATGCCAG",
         //||||||||||||||||||||||||||||||||||||          | |  |        ||      | |   | |   |    |  |  | |    | | || |        |  |        |    |  |  |    | |
          "AGATAAGTCCATGAAGTCACCAGCACCGTCCATGTTTCTCACTGCTTCCTCGGCGTTCCTCCAGAACCAAGCGTTACACCCCAACACAGGATGTGTGCCATAAATACTGGTTGCATGAATGGCTATTTTTTTTTAACTTCACTTTTTTCT",
         standardAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("36M114S"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(36U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getFStrandReferencePosition());
}

void TestSequencingAdapter::testStdBeforeSequence()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
        //CTGTCTCTTATACACATCT*
       //  ||||||||||||||||||-----------------------------------------------------------------------------------------------
    align("TGTCTCTTATACACATCTCCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAACACAACCATCGAGTCCACATCAGATATGCCAG",
         //          | |  |        ||      | |   | |   |    |  |  | |    | | || |        |  |        |    |  |  |    | |
          "CTCACTGCTTCCTCGGCGTTCCTCCAGAACCAAGCGTTACACCCCAACACAGGATGTGTGCCATAAATACTGGTTGCATGAATGGCTATTTTTTTTTAACTTCACTTTTTTCT",
         standardAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("113M"), fragmentMetadata.getCigarString());
}

void TestSequencingAdapter::testStdReverseAfterSequence()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = true;
                                             //                    *AGATGTGTATAAGAGACAG
       //  ---------------------------------------------------------||||||||||||||||||
    align("TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACA",
         //||||||x||x|||x|||x|||x|||x|||x||x|||x||xx|x|xx|x|xx|xxxxx|x|xxxxx|x||x|x|x|
          "TGGTTACGGCAGCTGTATAAGTGTGCTACTGCCATGCTCCCTTTTCTCTGGAATATGATAAAAAAAAAAAAAAAA",
         standardAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("75M"), fragmentMetadata.getCigarString());
}

void TestSequencingAdapter::testStdReverseSequenceTooGood()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = true;
                                             //                    *AGATGTGTATAAGAGACAG
       //  ---------------------------------------------------------|||||||||||||||||||
    align("TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG",
         //||||||x|||x|||x|||x|||x|||x||||||||||||| |||  | |  |     | |     | || | | |
          "TGGTTACGGTCGCGATAATAGCTTGTCACCGCTATGTTCTCTCTTCTCTGGAATATGATAAAAAAAAAAAAAAAAA",
         standardAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("76M"), fragmentMetadata.getCigarString());
}

/*
not supported cases:
original CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG

                                            CTGTCTCTTATACACATCT
                                                            AGATGTGTATAAGAGACAG
ATCATCACCTTAACCTCTGATAATCGTCATATACCGGACAAGACCTGTCTCTTATACACAAGATGTGTATAAGAGACAGGCTGGATACGTTGCAAAACAT

very short cluster with couple of bases missing in adapter:
                                    CTGTCTCTTATACACATCT
                                                     AGATGTGTATAAGAGACAG
AAGACTTTCACGCCTTCTTCAAACTCGGTCACTGGCCTGTCTCTTATACACATCTATGTGTATAAGAGACAGATCTGGGCCAGCTTCATGACATAACGCG
             CTTCTTCAAACTCGGTCACTGGCCTGTCTCTTATACACATCTATGTGTATAAGAGACAGATCTGGGCCAGCTTCATGACATAACGCGGTTGTTGAGATAA


clips off wrong side
                                CTGTCTCTTATACACATCTAGATGTGT AT   AAGAG ACA G
GTCGAATTGTGCGGGTAGCGATGCCATAAGCCCTGTCTCTTATACACATCTAGATGTGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTAT

One base mismatches towards the end of the adapter
                                         AGATGTGTATAAGAGACAG
<GACCACCGAGATCTACACTATCCTCTTCGTCGGCAGCGTCAGATGCGTATAAGAGACAGGCGCCAAACTTCGCCTACGAGTGGGCCGCACAGCGTGGACT

                                                CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG
CATGAATATTGTACGGTACCATAAATACTTGACCACCTGTAGTACATGAAAACCCAATCCACATCAAAACCCCCTCCCCATGCTTACAAGCAAGACCAGC>


unknown-flowcell_0:1:1:260110:0 147 chr1    7601960 0   114S36M =   7602015 111 CTGTCCTGATACACTGTGATGCTGTTTTTTTTTTTTTTTTTTTTAATGATCCGGCGACCACCGAGATCTACACAAGGAGTATCGTCGGCAGCGTCAGATGTGTATAAGAGACAGAGATAAGTCCATGAAGTCACCAGCACCGTCCATGTT  ((+(((+(+(+(+((((+((+++&&&&&)&0&)&.B>5BA:4+(:+20&&&<<.0300(2<<@@:3;@@9DCDDCCC;DDDDFFHEJIJJJJJJJJJJIJJJIJJJJJJJIJIJIJJJJJJJJJJJJJJIGHEJJJIHFHHHFFFFFCB@  SM:i:0  AS:i:0  RG:Z:0  NM:i:0  BC:Z:none
unknown-flowcell_0:1:1:260110:0 99  chr1    7602015 0   55S57M38S   =   7601960 -111    AGATAAGTCCATGAAGTCACCAGCACCGTCCATGTTCTGTCTCTTATACACATCTCCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAACACAACCATCGAGTCCACATCAGATATGCCAG  CCCFFFFFHHHHHJJJGHJJJJJJJJJJIJJJJJIJJJJJJJJJJJJJJJJJJIJJJJJJJJJJJHHFFDDDDDDDDDDDDDDDBDDDDDDDDBDDDCDDDCACDDDD>BD.0&))&)&+(+(((+(((((+(((+((+44(((((((((  SM:i:0  AS:i:0  RG:Z:0  NM:i:42 BC:Z:none

                                    CTGTCTCTTATACACATCT
                                    |||||||||||||||||||
AGATAAGTCCATGAAGTCACCAGCACCGTCCATGTTCTGTCTCTTATACACATCTCCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAACACAACCATCGAGTCCACATCAGATATGCCAG
AGATAAGTCCATGAAGTCACCAGCACCGTCCATGTTTCTCACTGCTTCCTCGGCGTTCCTCCAGAACCAAGCGTTACACCCCAACACAGGATGTGTGCCATAAATACTGGTTGCATGAATGGCTATTTTTTTTTAACTTCACTTTTTTCTTAATTAAAAATTT

*/



