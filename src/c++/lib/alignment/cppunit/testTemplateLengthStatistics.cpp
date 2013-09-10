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
#include "testTemplateLengthStatistics.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestTemplateLengthStatistics, registryName("TemplateLengthStatistics"));


TestTemplateLengthStatistics::TestTemplateLengthStatistics()
{
}

void TestTemplateLengthStatistics::setUp()
{
}

void TestTemplateLengthStatistics::tearDown()
{
}

void TestTemplateLengthStatistics::testAlignmentModels()
{
    using isaac::alignment::TemplateLengthStatistics;
    using isaac::alignment::FragmentMetadata;
    const std::string &(*alignmentModelName)(TemplateLengthStatistics::AlignmentModel) = &TemplateLengthStatistics::alignmentModelName;
    TemplateLengthStatistics::AlignmentModel (*alignmentModel)(const FragmentMetadata &, const FragmentMetadata &) = &TemplateLengthStatistics::alignmentModel;
    CPPUNIT_ASSERT_EQUAL(std::string("unknown"), alignmentModelName(TemplateLengthStatistics::InvalidAlignmentModel));
    isaac::alignment::FragmentMetadata f1, f2;
    f1.position = 0;
    f1.reverse = false;
    f2.position = 1;
    f2.reverse = false;
    CPPUNIT_ASSERT_EQUAL(std::string("FF+"), alignmentModelName(alignmentModel(f1, f2)));
    f2.reverse = true;
    CPPUNIT_ASSERT_EQUAL(std::string("FR+"), alignmentModelName(alignmentModel(f1, f2)));
    f1.reverse = true;
    CPPUNIT_ASSERT_EQUAL(std::string("RR+"), alignmentModelName(alignmentModel(f1, f2)));
    f2.reverse = false;
    CPPUNIT_ASSERT_EQUAL(std::string("RF+"), alignmentModelName(alignmentModel(f1, f2)));
    f1.position = 2;
    f1.reverse = false;
    f2.position = 1;
    f2.reverse = false;
    CPPUNIT_ASSERT_EQUAL(std::string("FF-"), alignmentModelName(alignmentModel(f1, f2)));
    f2.reverse = true;
    CPPUNIT_ASSERT_EQUAL(std::string("FR-"), alignmentModelName(alignmentModel(f1, f2)));
    f1.reverse = true;
    CPPUNIT_ASSERT_EQUAL(std::string("RR-"), alignmentModelName(alignmentModel(f1, f2)));
    f2.reverse = false;
    CPPUNIT_ASSERT_EQUAL(std::string("RF-"), alignmentModelName(alignmentModel(f1, f2)));
}

void TestTemplateLengthStatistics::testAlignmentClassNames()
{
    using isaac::alignment::TemplateLengthStatistics;
    using isaac::alignment::FragmentMetadata;
    TemplateLengthStatistics::AlignmentClass (*alignmentClass)(TemplateLengthStatistics::AlignmentModel) = &TemplateLengthStatistics::alignmentClass;
    const std::string &(*alignmentClassName)(TemplateLengthStatistics::AlignmentClass) = &TemplateLengthStatistics::alignmentClassName;
    CPPUNIT_ASSERT_EQUAL(std::string("F+"), alignmentClassName(alignmentClass(TemplateLengthStatistics::FFp)));
    CPPUNIT_ASSERT_EQUAL(std::string("R+"), alignmentClassName(alignmentClass(TemplateLengthStatistics::FRp)));
    CPPUNIT_ASSERT_EQUAL(std::string("R-"), alignmentClassName(alignmentClass(TemplateLengthStatistics::RFp)));
    CPPUNIT_ASSERT_EQUAL(std::string("F-"), alignmentClassName(alignmentClass(TemplateLengthStatistics::RRp)));
    CPPUNIT_ASSERT_EQUAL(std::string("F-"), alignmentClassName(alignmentClass(TemplateLengthStatistics::FFm)));
    CPPUNIT_ASSERT_EQUAL(std::string("R-"), alignmentClassName(alignmentClass(TemplateLengthStatistics::FRm)));
    CPPUNIT_ASSERT_EQUAL(std::string("R+"), alignmentClassName(alignmentClass(TemplateLengthStatistics::RFm)));
    CPPUNIT_ASSERT_EQUAL(std::string("F+"), alignmentClassName(alignmentClass(TemplateLengthStatistics::RRm)));
    CPPUNIT_ASSERT_EQUAL(std::string("unknown"), alignmentClassName(TemplateLengthStatistics::InvalidAlignmentClass));
}
void TestTemplateLengthStatistics::addTemplates(isaac::alignment::TemplateLengthDistribution & tls)
{
    using isaac::alignment::FragmentMetadata;

    const std::vector<unsigned > cigarBuffer(1, 16);
    std::vector<std::vector<FragmentMetadata> > f(2, std::vector<FragmentMetadata>(1));
    f[0][0].contigId = 0;
    f[0][0].position = 0;
    f[0][0].observedLength = 1;
    f[0][0].reverse = false;
    f[0][0].cigarBuffer = &cigarBuffer;
    f[0][0].cigarOffset = 0;
    f[0][0].cigarLength = 1;
    f[1][0] = f[0][0];
    f[1][0].reverse = true;
    for (unsigned int i = 1; i < 10000; ++i)
    {
        CPPUNIT_ASSERT_EQUAL(false, tls.addTemplate(f));
        ++f[1][0].position;
    }
    // swap read 1 and read 2 to have the reverse model represented as well
    std::swap(f[0], f[1]);
    CPPUNIT_ASSERT_EQUAL(false, tls.addTemplate(f));
    CPPUNIT_ASSERT_EQUAL(14U, tls.getStatistics().getMin());
    CPPUNIT_ASSERT_EQUAL(5001U, tls.getStatistics().getMedian());
    CPPUNIT_ASSERT_EQUAL(9987U, tls.getStatistics().getMax());
    CPPUNIT_ASSERT_EQUAL(3414U, tls.getStatistics().getLowStdDev());
    CPPUNIT_ASSERT_EQUAL(3413U, tls.getStatistics().getHighStdDev());
    // swap read 1 and read 2 in the original configuration
    std::swap(f[0], f[1]);
    f[1][0].position = f[0][0].position;
    for (unsigned int i = 1; i < 10000; ++i)
    {
        CPPUNIT_ASSERT_EQUAL(false, tls.addTemplate(f));
        ++f[1][0].position;
    }
    CPPUNIT_ASSERT_EQUAL(true, tls.addTemplate(f));
}

void TestTemplateLengthStatistics::testStatistics()
{
    using isaac::alignment::TemplateLengthDistribution;
    using isaac::alignment::FragmentMetadata;
    TemplateLengthDistribution tls(-1);
    addTemplates(tls);
    CPPUNIT_ASSERT_EQUAL(14U, tls.getStatistics().getMin());
    CPPUNIT_ASSERT_EQUAL(5001U, tls.getStatistics().getMedian());
    CPPUNIT_ASSERT_EQUAL(9987U, tls.getStatistics().getMax());
    CPPUNIT_ASSERT_EQUAL(3414U, tls.getStatistics().getLowStdDev());
    CPPUNIT_ASSERT_EQUAL(3413U, tls.getStatistics().getHighStdDev());
}

void TestTemplateLengthStatistics::testMateDriftRange()
{
    using isaac::alignment::TemplateLengthDistribution;
    using isaac::alignment::FragmentMetadata;
    TemplateLengthDistribution tls(123);
    addTemplates(tls);
    CPPUNIT_ASSERT_EQUAL(5001U, tls.getStatistics().getMedian());

    CPPUNIT_ASSERT_EQUAL(tls.getStatistics().getMedian() - 123U, tls.getStatistics().getMateMin());
    CPPUNIT_ASSERT_EQUAL(tls.getStatistics().getMedian() + 123U, tls.getStatistics().getMateMax());
}


void TestTemplateLengthStatistics::testNoMateDriftRange()
{
    using isaac::alignment::TemplateLengthDistribution;
    using isaac::alignment::FragmentMetadata;
    TemplateLengthDistribution tls(-1);
    addTemplates(tls);
    CPPUNIT_ASSERT_EQUAL(5001U, tls.getStatistics().getMedian());
    CPPUNIT_ASSERT_EQUAL(tls.getStatistics().getMin(), tls.getStatistics().getMateMin());
    CPPUNIT_ASSERT_EQUAL(tls.getStatistics().getMax(), tls.getStatistics().getMateMax());
}


void TestTemplateLengthStatistics::testMateOrientation()
{
    using isaac::alignment::TemplateLengthStatistics;
    using isaac::alignment::FragmentMetadata;
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::FRp, TemplateLengthStatistics::RFm, -1);
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(0, false));
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(1, false));
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(0, true));
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(1, true));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::RFp, TemplateLengthStatistics::FRm, -1);
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(0, false));
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(1, false));
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(0, true));
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(1, true));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::FRm, TemplateLengthStatistics::RFp, -1);
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(0, false));
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(1, false));
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(0, true));
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(1, true));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::RFm, TemplateLengthStatistics::FRp, -1);
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(0, false));
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(1, false));
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(0, true));
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(1, true));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::FFm, TemplateLengthStatistics::RRp, -1);
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(0, false));
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(1, false));
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(0, true));
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(1, true));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::RRm, TemplateLengthStatistics::FFp, -1);
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(0, false));
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(1, false));
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(0, true));
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(1, true));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::FFp, TemplateLengthStatistics::RRm, -1);
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(0, false));
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(1, false));
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(0, true));
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(1, true));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::RRp, TemplateLengthStatistics::FFm, -1);
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(0, false));
        CPPUNIT_ASSERT_EQUAL(false, tls.mateOrientation(1, false));
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(0, true));
        CPPUNIT_ASSERT_EQUAL(true, tls.mateOrientation(1, true));
    }
}
void TestTemplateLengthStatistics::testMateMinPosition()
{
    using isaac::alignment::TemplateLengthStatistics;
    using isaac::alignment::FragmentMetadata;
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::FRp, TemplateLengthStatistics::RFm, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(517L, tls.mateMinPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(367L, tls.mateMinPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(533L, tls.mateMinPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(383L, tls.mateMinPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::RFp, TemplateLengthStatistics::FRm, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(367L, tls.mateMinPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(517L, tls.mateMinPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(383L, tls.mateMinPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(533L, tls.mateMinPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::FRm, TemplateLengthStatistics::RFp, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(367L, tls.mateMinPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(517L, tls.mateMinPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(383L, tls.mateMinPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(533L, tls.mateMinPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::RFm, TemplateLengthStatistics::FRp, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(517L, tls.mateMinPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(367L, tls.mateMinPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(533L, tls.mateMinPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(383L, tls.mateMinPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::FFp, TemplateLengthStatistics::RRm, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(517L, tls.mateMinPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(367L, tls.mateMinPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(383L, tls.mateMinPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(533L, tls.mateMinPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::RRp, TemplateLengthStatistics::FFm, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(367L, tls.mateMinPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(517L, tls.mateMinPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(533L, tls.mateMinPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(383L, tls.mateMinPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::FFm, TemplateLengthStatistics::RRp, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(367L, tls.mateMinPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(517L, tls.mateMinPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(533L, tls.mateMinPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(383L, tls.mateMinPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::RRm, TemplateLengthStatistics::FFp, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(517L, tls.mateMinPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(367L, tls.mateMinPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(383L, tls.mateMinPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(533L, tls.mateMinPosition(1, true, 500, readLengths));
    }
}

void TestTemplateLengthStatistics::testMateMaxPosition()
{
    using isaac::alignment::TemplateLengthStatistics;
    using isaac::alignment::FragmentMetadata;
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::FRp, TemplateLengthStatistics::RFm, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(617L, tls.mateMaxPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(467L, tls.mateMaxPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(633L, tls.mateMaxPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(483L, tls.mateMaxPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::RFp, TemplateLengthStatistics::FRm, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(467L, tls.mateMaxPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(617L, tls.mateMaxPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(483L, tls.mateMaxPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(633L, tls.mateMaxPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::FRm, TemplateLengthStatistics::RFp, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(467L, tls.mateMaxPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(617L, tls.mateMaxPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(483L, tls.mateMaxPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(633L, tls.mateMaxPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::RFm, TemplateLengthStatistics::FRp, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(617L, tls.mateMaxPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(467L, tls.mateMaxPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(633L, tls.mateMaxPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(483L, tls.mateMaxPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::FFp, TemplateLengthStatistics::RRm, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(617L, tls.mateMaxPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(467L, tls.mateMaxPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(483L, tls.mateMaxPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(633L, tls.mateMaxPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::RRp, TemplateLengthStatistics::FFm, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(467L, tls.mateMaxPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(617L, tls.mateMaxPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(633L, tls.mateMaxPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(483L, tls.mateMaxPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::FFm, TemplateLengthStatistics::RRp, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(467L, tls.mateMaxPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(617L, tls.mateMaxPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(633L, tls.mateMaxPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(483L, tls.mateMaxPosition(1, true, 500, readLengths));
    }
    {
        TemplateLengthStatistics tls(100, 200, 170, 160, 175, TemplateLengthStatistics::RRm, TemplateLengthStatistics::FFp, -1);
        const unsigned readLengths[] = {67, 83};
        CPPUNIT_ASSERT_EQUAL(617L, tls.mateMaxPosition(0, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(467L, tls.mateMaxPosition(0, true, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(483L, tls.mateMaxPosition(1, false, 500, readLengths));
        CPPUNIT_ASSERT_EQUAL(633L, tls.mateMaxPosition(1, true, 500, readLengths));
    }
}
