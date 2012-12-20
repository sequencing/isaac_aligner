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
 **
 ** \file testSemialignedClipper.cpp
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
#include "testSemialignedClipper.hh"

#include "alignment/Cluster.hh"
#include "alignment/matchSelector/FragmentSequencingAdapterClipper.hh"
#include "flowcell/SequencingAdapterMetadata.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestSemialignedClipper, registryName("SemialignedClipper"));

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

TestSemialignedClipper::TestSemialignedClipper() :
    readMetadataList(getReadMetadataList()),
    seedMetadataList(getSeedMetadataList()),
    flowcells(1, isaac::flowcell::Layout("", isaac::flowcell::Layout::Fastq, std::vector<unsigned>(),
                                         readMetadataList, seedMetadataList, "blah")),
    fragmentBuilder_(flowcells, 123, seedMetadataList.size()/2, 8, true)
{

}

void TestSemialignedClipper::setUp()
{
}

void TestSemialignedClipper::tearDown()
{
}

inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

namespace testSemialignedClipper
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

} // namespace testSemialignedClipper

namespace isaac
{
namespace alignment
{

inline void phredToBcl(char &qual)
{
    qual -= 33;
}

template<class InpuT> InpuT& operator >>(InpuT &input, isaac::alignment::Read &read);

template<> testSemialignedClipper::ReadInit& operator >><testSemialignedClipper::ReadInit >(
    testSemialignedClipper::ReadInit &input,
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

void TestSemialignedClipper::testEverything()
{
    testLeftClipForward();
    testRightClipForward();
    testRightClipStartBeforeRef();
    testLeftClipStartBeforeRef();

}

static const isaac::reference::Contig makeContig(const std::string forward, long &firstPosOffset)
{
    std::string::const_iterator begin = std::find_if(forward.begin(), forward.end(),
                                                     boost::bind(&boost::cref<char>, _1) != ' ');
    isaac::reference::Contig ret(0, "vasja");
    ret.forward_ = std::vector<char>(begin, forward.end());
    firstPosOffset = -std::distance(forward.begin(), begin);
    return ret;
}

void TestSemialignedClipper::align(
    const std::string &read,
    const std::string &reference,
    const isaac::alignment::matchSelector::SequencingAdapterList &adapters,
    isaac::alignment::FragmentMetadata &fragmentMetadata)
{
    isaac::alignment::Cluster cluster(isaac::flowcell::getMaxReadLength(flowcells));
    testSemialignedClipper::ReadInit init(read, fragmentMetadata.reverse);
    init >> cluster.at(0);

    if (fragmentMetadata.isNoMatch())
    {
        fragmentMetadata.contigId = 0;
        fragmentMetadata.position = 0;
    }
    fragmentMetadata.cluster = &cluster;
    fragmentMetadata.cigarBuffer = &cigarBuffer_;
    std::vector<isaac::reference::Contig> contigList;
    contigList.push_back(makeContig(reference, fragmentMetadata.position));

    isaac::alignment::matchSelector::FragmentSequencingAdapterClipper adapterClipper(adapters);
    adapterClipper.checkInitStrand(fragmentMetadata, contigList.at(0));
    fragmentBuilder_.alignUngapped(fragmentMetadata, cigarBuffer_, readMetadataList, adapterClipper, contigList.at(0));
    clipper_.clip(contigList, fragmentMetadata);
}


void TestSemialignedClipper::testLeftClipForward()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    align("AGATCTACACATATCCGCCACGTGGACAGAGAATATGTGTAGATCTACACATATTCTCTGTCTTGTAACGCCATTGTGCGAAAATGGCGATGGAATTGGT",
         //|x|xxx|x|x|x|x||||||||||x|x|x|x||x|xxxxx|x|xxxxx|xx|x|xx|||x||||||||||||||||||||||||||||||||||||||||
          "AAAAAAAAAAAAAACCGCCACGTGAAAAAAAAAAAAAAAAATAGTATAAGGTCTCTTCTCTCTTGTAACGCCATTGTGCGAAAATGGCGATGGAATTGGT",
         noAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("14S86M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 14U), fragmentMetadata.getStrandReferencePosition());
}

void TestSemialignedClipper::testRightClipForward()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    align("TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         //||||||||||||||||||||||||||||||||||||||||x|||xx|x|xx|xxxxx|x|xxxxx|x||x|x|x|x||||||||||x|x|x|x|xxx|x|
          "TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTCTCTTCTCTGGAATATGATAAAAAAAAAAAAAAAAAGTGCACCGCCAAAAAAAAAAAAAA",
         noAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("86M14S"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getStrandReferencePosition());
}

void TestSemialignedClipper::testRightClipStartBeforeRef()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    align("AAAAATGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACAT",
         //     ||||||||||||||||||||||||||||||||||||||||x|||xx|x|xx|xxxxx|x|xxxxx|x||x|x|x|x||||||||||x|x|x|x|xxx|x|
          "     TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTCTCTTCTCTGGAATATGATAAAAAAAAAAAAAAAAAGTGCACCGCCAAAAAAAAAAAAAA",
         noAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getStrandReferencePosition());
    CPPUNIT_ASSERT_EQUAL(std::string("5S86M9S"), fragmentMetadata.getCigarString());
}

void TestSemialignedClipper::testLeftClipStartBeforeRef()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    align("AAAAAGATCTACACATATCCGCCACGTGGACAGAGAATATGTGTAGATCTACACATATTCTCTGTCTTGTAACGCCATTGTGCGAAAATGGCGATGGAAT",
         //    |x|xxx|x|x|x|x||||||||||x|x|x|x||x|xxxxx|x|xxxxx|xx|x|xx|||x||||||||||||||||||||||||||||||||||||||||
          "    AAAAAAAAAAAAAACCGCCACGTGAAAAAAAAAAAAAAAAATAGTATAAGGTCTCTTCTCTCTTGTAACGCCATTGTGCGAAAATGGCGATGGAATTGGT",
         noAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("18S82M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 14U), fragmentMetadata.getStrandReferencePosition());
}


