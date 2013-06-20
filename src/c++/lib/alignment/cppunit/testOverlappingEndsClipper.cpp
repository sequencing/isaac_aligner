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
 ** \file testOverlappingEndsClipper.cpp
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
#include "testOverlappingEndsClipper.hh"

#include "alignment/Cluster.hh"
#include "alignment/matchSelector/OverlappingEndsClipper.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestOverlappingEndsClipper, registryName("OverlappingEndsClipper"));

TestOverlappingEndsClipper::TestOverlappingEndsClipper() : cluster_(1234)
{

}

void TestOverlappingEndsClipper::setUp()
{
}

void TestOverlappingEndsClipper::tearDown()
{
}

inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

namespace testOverlappingEndsClipper
{

//static const std::string irrelevantQualities("CFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDB<F>FGGBFGGFGCGGGDGGDDFHHHFEGGBGDGGBGGBEGEGGBGEHDHHHGGGGGDGGGG?GGGG");

struct ReadInit : public std::pair<std::string, std::string>
{
    static std::string reverseString(std::string fwd)
    {
        std::reverse(fwd.begin(), fwd.end());
        return fwd;
    }

    typedef std::pair<std::string, std::string> BaseType;
    ReadInit(const std::string &read, const std::string &quality, const bool reverse) :
        BaseType(reverse ? reverseString(read) : read, reverse ? quality : reverseString(quality))
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

template<> testOverlappingEndsClipper::ReadInit& operator >><testOverlappingEndsClipper::ReadInit >(
    testOverlappingEndsClipper::ReadInit &input,
    isaac::alignment::Read &read)
{
    ISAAC_THREAD_CERR << input.first << " " << input.second << std::endl;
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

void TestOverlappingEndsClipper::testEverything()
{
    isaac::reference::ContigList contigList;
    isaac::alignment::BamTemplate templ(cigarBuffer_);


    {
        init("ACGT", "CFCE", false,
             " ACGT", "BDBE", true,
             "ACGT", templ, contigList);

        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(1L, templ.getFragmentMetadata(1).position);

        isaac::alignment::matchSelector::OverlappingEndsClipper clipper;
        clipper.clip(contigList, templ);

        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("3S1M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(4L, templ.getFragmentMetadata(1).position);
    }

    {
        init("ACGT", "BAAA", false,
             " ACGT", "CFCE", true,
             "ACGT", templ, contigList);

        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(1L, templ.getFragmentMetadata(1).position);

        isaac::alignment::matchSelector::OverlappingEndsClipper clipper;
        clipper.clip(contigList, templ);

        CPPUNIT_ASSERT_EQUAL(std::string("1M3S"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(1L, templ.getFragmentMetadata(1).position);
    }


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

inline isaac::flowcell::ReadMetadataList getReadMetadataList(const unsigned l0, const unsigned l1)
{
    std::vector<isaac::flowcell::ReadMetadata> ret =
        boost::assign::list_of
            (isaac::flowcell::ReadMetadata(1, l0, 0, 0))
            (isaac::flowcell::ReadMetadata(l0 + 1, l0 + l1, 1, l0))
            ;
    return ret;
}

void TestOverlappingEndsClipper::init(
    const std::string &read1, const std::string &quality1, const bool read1Reverse,
    const std::string &read2, const std::string &quality2, const bool read2Reverse,
    const std::string &reference,
    isaac::alignment::BamTemplate &templ,
    isaac::reference::ContigList &contigList)
{
    contigList.clear();
    long firstPosOffset = 0;
    contigList.push_back(makeContig(reference, firstPosOffset));

    const unsigned r1Start = read1.find_first_not_of(' ');
    const unsigned r2Start = read2.find_first_not_of(' ');
    static isaac::flowcell::ReadMetadataList readMetadataList;
    readMetadataList = getReadMetadataList(read1.length() - r1Start, read2.length() - r2Start);

    testOverlappingEndsClipper::ReadInit init1(read1.substr(r1Start), quality1, read1Reverse);
    init1 >> cluster_.at(0);
    testOverlappingEndsClipper::ReadInit init2(read2.substr(r2Start), quality2, read2Reverse);
    init2 >> cluster_.at(1);

    templ.initialize(readMetadataList, cluster_);

    cigarBuffer_.clear();

    templ.getFragmentMetadata(0).reverse = read1Reverse;
    templ.getFragmentMetadata(0).cigarBuffer = &cigarBuffer_;
    templ.getFragmentMetadata(0).cigarOffset = cigarBuffer_.size();
    cigarBuffer_.push_back(isaac::alignment::Cigar::encode(readMetadataList.at(0).getLength(), isaac::alignment::Cigar::ALIGN));
    templ.getFragmentMetadata(0).cigarLength = cigarBuffer_.size() - templ.getFragmentMetadata(0).cigarOffset;
    templ.getFragmentMetadata(0).contigId = 0;
    templ.getFragmentMetadata(0).position = read1.find_first_not_of(' ');
    templ.getFragmentMetadata(0).observedLength = readMetadataList.at(0).getLength();

    templ.getFragmentMetadata(1).reverse = read2Reverse;
    templ.getFragmentMetadata(1).cigarBuffer = &cigarBuffer_;
    templ.getFragmentMetadata(1).cigarOffset = cigarBuffer_.size();
    cigarBuffer_.push_back(isaac::alignment::Cigar::encode(readMetadataList.at(1).getLength(), isaac::alignment::Cigar::ALIGN));
    templ.getFragmentMetadata(1).cigarLength = cigarBuffer_.size() - templ.getFragmentMetadata(1).cigarOffset;
    templ.getFragmentMetadata(1).contigId = 0;
    templ.getFragmentMetadata(1).position = read2.find_first_not_of(' ');
    templ.getFragmentMetadata(1).observedLength = readMetadataList.at(1).getLength();
}


