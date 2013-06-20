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
 ** \file testFragmentBuilder2.cpp
 **
 ** More fragment builder tests.
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
#include "testFragmentBuilder2.hh"

#include "alignment/fragmentBuilder/GappedAligner.hh"
#include "alignment/fragmentBuilder/UngappedAligner.hh"
#include "alignment/matchSelector/FragmentSequencingAdapterClipper.hh"
#include "alignment/BandedSmithWaterman.hh"
#include "alignment/Cluster.hh"
#include "flowcell/SequencingAdapterMetadata.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestFragmentBuilder2, registryName("FragmentBuilder2"));

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

static const isaac::alignment::matchSelector::SequencingAdapterList noAdapters;

TestFragmentBuilder2::TestFragmentBuilder2() :
    readMetadataList(getReadMetadataList()),
    seedMetadataList(getSeedMetadataList()),
    flowcells(1, isaac::flowcell::Layout("", isaac::flowcell::Layout::Fastq, std::vector<unsigned>(),
                                         readMetadataList, seedMetadataList, "blah"))
{

}

void TestFragmentBuilder2::setUp()
{
}

void TestFragmentBuilder2::tearDown()
{
}

inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

namespace testFragmentBuilder2
{

static const std::string irrelevantQualities("CFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDB<F>FGGBFGGFGCGGGDGGDDFHHHFEGGBGDGGBGGBEGEGGBGEHDHHHGGGGGDGGGG?GGGG");

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

    }
};

} // namespace testFragmentBuilder2

namespace isaac
{
namespace alignment
{

inline void phredToBcl(char &qual)
{
    qual -= 33;
}

template<class InpuT> InpuT& operator >>(InpuT &input, isaac::alignment::Read &read);

template<> testFragmentBuilder2::ReadInit& operator >><testFragmentBuilder2::ReadInit >(
    testFragmentBuilder2::ReadInit &input,
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

void TestFragmentBuilder2::testEverything()
{
    testMismatchCount();
    testMismatchCycles();
    testMismatchCyclesWithSoftClip();
    testGapped();
    testGappedWithNs();
}

static const isaac::reference::Contig makeContig(const std::string forward)
{
    isaac::reference::Contig ret(0, "vasja");
    ret.forward_ = vectorFromString(forward);
    return ret;
}

static const int ELAND_MATCH_SCORE = 2;
static const int ELAND_MISMATCH_SCORE = -1;
static const int ELAND_GAP_OPEN_SCORE = -15;
static const int ELAND_GAP_EXTEND_SCORE = -3;
static const int ELAND_MIN_GAP_EXTEND_SCORE = 25;

void TestFragmentBuilder2::align(
    const std::string &read,
    const std::string &reference,
    const isaac::alignment::matchSelector::SequencingAdapterList &adapters,
    isaac::alignment::FragmentMetadata &fragmentMetadata,
    const bool gapped)
{
    isaac::alignment::Cluster cluster(isaac::flowcell::getMaxReadLength(flowcells));
    testFragmentBuilder2::ReadInit init(read, fragmentMetadata.reverse);
    init >> cluster.at(0);

    if (fragmentMetadata.isNoMatch())
    {
        fragmentMetadata.contigId = 0;
        fragmentMetadata.position = 0;
    }
    fragmentMetadata.cluster = &cluster;
    fragmentMetadata.cigarBuffer = &cigarBuffer_;

    const isaac::reference::Contig referenceContig = makeContig(reference);

    isaac::alignment::matchSelector::FragmentSequencingAdapterClipper adapterClipper(adapters);
    adapterClipper.checkInitStrand(fragmentMetadata, referenceContig);

    isaac::alignment::fragmentBuilder::UngappedAligner ungappedAligner(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE);

    ungappedAligner.alignUngapped(fragmentMetadata, cigarBuffer_, readMetadataList, adapterClipper, referenceContig);
    if (gapped)
    {
        isaac::alignment::fragmentBuilder::GappedAligner gappedAligner(flowcells, false, ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE);
        isaac::alignment::FragmentMetadata tmp = fragmentMetadata;
        const unsigned matchCount = gappedAligner.alignGapped(tmp, cigarBuffer_, readMetadataList, adapterClipper, referenceContig);
        if (matchCount + isaac::alignment::BandedSmithWaterman::WIDEST_GAP_SIZE > fragmentMetadata.getObservedLength() &&
                                (tmp.mismatchCount <= 5) &&
                                (fragmentMetadata.mismatchCount > tmp.mismatchCount) &&
                                fragmentMetadata.logProbability < tmp.logProbability)
        {
            fragmentMetadata = tmp;
        }
    }
}


void TestFragmentBuilder2::testMismatchCount()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    align("TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         //||||||||||||||||||||||||||||||||||||||||x|||xx|x|xx|xxxxx|x|xxxxx|x||x|x|x|x||||||||||x|x|x|x|xxx|x|
          "TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTCTCTTCTCTGGAATATGATAAAAAAAAAAAAAAAAAGTGCACCGCCAAAAAAAAAAAAAA",
         noAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("100M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(30U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(30U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(100U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getStrandReferencePosition());
}

void TestFragmentBuilder2::testMismatchCycles()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = true;

    align("TGGTTAAGATAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         //||||||||x|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          "TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         noAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("100M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(100U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getFStrandReferencePosition());
    CPPUNIT_ASSERT_EQUAL(92U, unsigned(*fragmentMetadata.getMismatchCyclesBegin()));
}

void TestFragmentBuilder2::testMismatchCyclesWithSoftClip()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    fragmentMetadata.contigId = 0;
    fragmentMetadata.position = -2;

    align("TTTGGTTAAGATAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTA",
         //**||||||||x|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            "TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         noAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("2S98M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(98U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getFStrandReferencePosition());
    CPPUNIT_ASSERT_EQUAL(11U, unsigned(*fragmentMetadata.getMismatchCyclesBegin()));
}

void TestFragmentBuilder2::testGapped()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    fragmentMetadata.contigId = 0;
    fragmentMetadata.position = 1;

    align( "TTTGGTTAAGATAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTA",
         // ||||||||x||||||||||||||||||||||||||||||||||||||||||||^\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
          "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         noAdapters,
         fragmentMetadata, true);

    CPPUNIT_ASSERT_EQUAL(std::string("53M1D47M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(2U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(101U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 1U), fragmentMetadata.getFStrandReferencePosition());
}

void TestFragmentBuilder2::testGappedWithNs()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    fragmentMetadata.contigId = 0;
    fragmentMetadata.position = 1;

    align( "TTTGGTTAAGATAGCGGTAAAAGCGTGTTACCGCAATGTTCTGnnnnTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTA",
         // |||||||||||||||||||||||||||||||||||||||||||||||||||||^\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
         // ||||||||x|||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||
          "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         noAdapters,
         fragmentMetadata, true);

    CPPUNIT_ASSERT_EQUAL(std::string("53M1D47M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(6U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(101U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 1U), fragmentMetadata.getFStrandReferencePosition());
}
