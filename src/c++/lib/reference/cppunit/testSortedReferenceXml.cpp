/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** BSD 2-Clause License
 **
 ** You should have received a copy of the BSD 2-Clause License
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testSortedReferenceXml.hh"

#include "xml/XmlReader.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestSortedReferenceXml, registryName("SortedReferenceXml"));

TestSortedReferenceXml::TestSortedReferenceXml()
    : xmlString(
"<?xml version=\"1.0\"?>\n"
"<SortedReference>\n"
"  <FormatVersion>2</FormatVersion>\n"
"  <SoftwareVersion>iSAAC-01.12.12.12</SoftwareVersion>"
"  <Contigs>\n"
"    <Contig Position=\"0\">\n"
"      <Index>0</Index>\n"
"      <KaryotypeIndex>1</KaryotypeIndex>"
"      <Name>10</Name>\n"
"      <Sequence>\n"
"        <File>/illumina/scratch/smallprojects/iSAAC/Genomes/HumanNCBI36Fasta_all.fa</File>\n"
"        <Offset>54</Offset>\n"
"        <Size>137630983</Size>\n"
"      </Sequence>\n"
"      <TotalBases>135374737</TotalBases>\n"
"      <AcgtBases>131624728</AcgtBases>\n"
"      <BamMetadata>"
"        <Sq>"
"          <As>Tada</As>"
"          <Ur>/blah</Ur>"
"          <M5>12345</M5>"
"        </Sq>"
"      </BamMetadata>"
"    </Contig>\n"
"    <Contig Position=\"135374737\">\n"
"      <Index>1</Index>\n"
"      <Name>11</Name>\n"
"      <Sequence>\n"
"        <File>/illumina/scratch/smallprojects/iSAAC/Genomes/HumanNCBI36Fasta_all.fa</File>\n"
"        <Offset>137631091</Offset>\n"
"        <Size>136693258</Size>\n"
"      </Sequence>\n"
"      <TotalBases>134452384</TotalBases>\n"
"      <AcgtBases>131130753</AcgtBases>\n"
"      <BamMetadata>"
"        <Sq>"
"          <As/>"
"          <M5>123456</M5>"
"        </Sq>"
"      </BamMetadata>"
"    </Contig>\n"
"  </Contigs>\n"
"  <Permutations>\n"
"    <Permutation Name=\"ABCD\">\n"
"      <Masks Width=\"1\">\n"
"        <Mask Mask=\"0\">\n"
"          <File>/illumina/scratch/smallprojects/iSAAC/Genomes/HumanNCBI36-neighbors/HumanNCBI36Fasta_all.fa-32mer-6bit-ABCD-00.dat</File>\n"
"          <Kmers>\n"
"            <Total>107990454</Total>\n"
"          </Kmers>\n"
"          <MaxPrefixRangeCount>927297</MaxPrefixRangeCount>\n"
"        </Mask>\n"
"        <Mask Mask=\"1\">\n"
"          <File>/illumina/scratch/smallprojects/iSAAC/Genomes/HumanNCBI36-neighbors/HumanNCBI36Fasta_all.fa-32mer-6bit-ABCD-01.dat</File>\n"
"          <Kmers>\n"
"            <Total>40727835</Total>\n"
"          </Kmers>\n"
"          <MaxPrefixRangeCount>42468</MaxPrefixRangeCount>\n"
"        </Mask>\n"
"        <Mask Mask=\"2\">\n"
"          <File>/illumina/scratch/smallprojects/iSAAC/Genomes/HumanNCBI36-neighbors/HumanNCBI36Fasta_all.fa-32mer-6bit-ABCD-02.dat</File>\n"
"          <Kmers>\n"
"            <Total>56150179</Total>\n"
"          </Kmers>\n"
"          <MaxPrefixRangeCount>103275</MaxPrefixRangeCount>\n"
"        </Mask>\n"
"      </Masks>\n"
"    </Permutation>\n"
"  </Permutations>\n"
"</SortedReference>\n"
)
{
}

void TestSortedReferenceXml::setUp()
{
}

void TestSortedReferenceXml::tearDown()
{
}

void TestSortedReferenceXml::checkContent(const isaac::reference::SortedReferenceMetadata &sortedReferenceMetadata)
{
    checkContigs(sortedReferenceMetadata);
    checkMasks(sortedReferenceMetadata);
}

void TestSortedReferenceXml::checkContigs(const isaac::reference::SortedReferenceMetadata &sortedReferenceMetadata)
{
    CPPUNIT_ASSERT_EQUAL(2U, unsigned(sortedReferenceMetadata.getContigsCount()));

    CPPUNIT_ASSERT_EQUAL(0UL, sortedReferenceMetadata.getContigs().at(0).genomicPosition_);
    CPPUNIT_ASSERT_EQUAL(0U, sortedReferenceMetadata.getContigs().at(0).index_);
    CPPUNIT_ASSERT_EQUAL(std::string("10"), sortedReferenceMetadata.getContigs().at(0).name_);
    CPPUNIT_ASSERT_EQUAL(boost::filesystem::path("/illumina/scratch/smallprojects/iSAAC/Genomes/HumanNCBI36Fasta_all.fa"), sortedReferenceMetadata.getContigs().at(0).filePath_);
    CPPUNIT_ASSERT_EQUAL(54UL, sortedReferenceMetadata.getContigs().at(0).offset_);
    CPPUNIT_ASSERT_EQUAL(137630983UL, sortedReferenceMetadata.getContigs().at(0).size_);
    CPPUNIT_ASSERT_EQUAL(135374737UL, sortedReferenceMetadata.getContigs().at(0).totalBases_);
    CPPUNIT_ASSERT_EQUAL(131624728UL, sortedReferenceMetadata.getContigs().at(0).acgtBases_);
    CPPUNIT_ASSERT_EQUAL(1U, sortedReferenceMetadata.getContigs().at(0).karyotypeIndex_);
    CPPUNIT_ASSERT_EQUAL(std::string("Tada"), sortedReferenceMetadata.getContigs().at(0).bamSqAs_);
    CPPUNIT_ASSERT_EQUAL(std::string("/blah"), sortedReferenceMetadata.getContigs().at(0).bamSqUr_);
    CPPUNIT_ASSERT_EQUAL(std::string("12345"), sortedReferenceMetadata.getContigs().at(0).bamM5_);

    CPPUNIT_ASSERT_EQUAL(135374737UL, sortedReferenceMetadata.getContigs().at(1).genomicPosition_);
    CPPUNIT_ASSERT_EQUAL(1U, sortedReferenceMetadata.getContigs().at(1).index_);
    CPPUNIT_ASSERT_EQUAL(std::string("11"), sortedReferenceMetadata.getContigs().at(1).name_);
    CPPUNIT_ASSERT_EQUAL(boost::filesystem::path("/illumina/scratch/smallprojects/iSAAC/Genomes/HumanNCBI36Fasta_all.fa"), sortedReferenceMetadata.getContigs().at(1).filePath_);
    CPPUNIT_ASSERT_EQUAL(137631091UL, sortedReferenceMetadata.getContigs().at(1).offset_);
    CPPUNIT_ASSERT_EQUAL(136693258UL, sortedReferenceMetadata.getContigs().at(1).size_);
    CPPUNIT_ASSERT_EQUAL(134452384UL, sortedReferenceMetadata.getContigs().at(1).totalBases_);
    CPPUNIT_ASSERT_EQUAL(131130753UL, sortedReferenceMetadata.getContigs().at(1).acgtBases_);
    CPPUNIT_ASSERT_EQUAL(1U, sortedReferenceMetadata.getContigs().at(1).karyotypeIndex_);
    CPPUNIT_ASSERT_EQUAL(std::string(""), sortedReferenceMetadata.getContigs().at(1).bamSqAs_);
    CPPUNIT_ASSERT_EQUAL(std::string(""), sortedReferenceMetadata.getContigs().at(1).bamSqUr_);
    CPPUNIT_ASSERT_EQUAL(std::string("123456"), sortedReferenceMetadata.getContigs().at(1).bamM5_);
}

void TestSortedReferenceXml::checkMasks(const isaac::reference::SortedReferenceMetadata &sortedReferenceMetadata)
{
    const isaac::reference::SortedReferenceMetadata::MaskFiles &list = sortedReferenceMetadata.getMaskFileList(32);
    CPPUNIT_ASSERT_EQUAL(3U, unsigned(list.size()));
    CPPUNIT_ASSERT_EQUAL(boost::filesystem::path("/illumina/scratch/smallprojects/iSAAC/Genomes/HumanNCBI36-neighbors/HumanNCBI36Fasta_all.fa-32mer-6bit-ABCD-02.dat"),
                         list.back().path);
    CPPUNIT_ASSERT_EQUAL(sortedReferenceMetadata.getDefaultMaskWidth(), list.back().maskWidth);
    CPPUNIT_ASSERT_EQUAL(1U, sortedReferenceMetadata.getDefaultMaskWidth());
}


void TestSortedReferenceXml::testAll()
{
    std::istringstream is(xmlString);
    isaac::reference::SortedReferenceMetadata sortedReferenceMetadata = isaac::reference::loadSortedReferenceXml(is);
    checkContent(sortedReferenceMetadata);
}

void TestSortedReferenceXml::testWriter()
{
    std::istringstream is(xmlString);
    isaac::reference::SortedReferenceMetadata sortedReferenceMetadata = isaac::reference::loadSortedReferenceXml(is);

    std::ostringstream os;
    isaac::reference::saveSortedReferenceXml(os, sortedReferenceMetadata);

//    ISAAC_THREAD_CERR << os.str() << std::endl;
    std::istringstream is2(os.str());
    checkContent(isaac::reference::loadSortedReferenceXml(is2));
}

void TestSortedReferenceXml::testContigsOnly()
{
    std::istringstream is(xmlString);
    isaac::reference::SortedReferenceMetadata sortedReferenceMetadata = isaac::reference::loadSortedReferenceXml(is);
    CPPUNIT_ASSERT_EQUAL(true, sortedReferenceMetadata.supportsSeedLength(32));
    sortedReferenceMetadata.clearMasks();

    std::ostringstream os;
    isaac::reference::saveSortedReferenceXml(os, sortedReferenceMetadata);

//    ISAAC_THREAD_CERR << os.str() << std::endl;

    std::istringstream is2(os.str());
    checkContigs(isaac::reference::loadSortedReferenceXml(is2));

}


void TestSortedReferenceXml::testMasksOnly()
{
    std::istringstream is(xmlString);
    isaac::reference::SortedReferenceMetadata sortedReferenceMetadata = isaac::reference::loadSortedReferenceXml(is);
    sortedReferenceMetadata.getContigs().clear();

    std::ostringstream os;
    isaac::reference::saveSortedReferenceXml(os, sortedReferenceMetadata);

//    ISAAC_THREAD_CERR << os.str() << std::endl;

    std::istringstream is2(os.str());
    checkMasks(isaac::reference::loadSortedReferenceXml(is2));
}

void TestSortedReferenceXml::testMerge()
{
    std::istringstream is1(xmlString);
    isaac::reference::SortedReferenceMetadata contigsReference = isaac::reference::loadSortedReferenceXml(is1);
    contigsReference.clearMasks();

    std::istringstream is2(xmlString);
    isaac::reference::SortedReferenceMetadata mask0Reference = isaac::reference::loadSortedReferenceXml(is2);
    mask0Reference.getContigs().clear();
    mask0Reference.getMaskFileList(32).pop_back();
    mask0Reference.getMaskFileList(32).pop_back();

    std::istringstream is3(xmlString);
    isaac::reference::SortedReferenceMetadata mask1Reference = isaac::reference::loadSortedReferenceXml(is3);
    mask1Reference.getContigs().clear();
    mask1Reference.getMaskFileList(32).erase(mask1Reference.getMaskFileList(32).begin());

    isaac::reference::SortedReferenceMetadata mergedReference;
    mergedReference.merge(contigsReference);
    mergedReference.merge(mask0Reference);
    mergedReference.merge(mask1Reference);

    checkContent(mergedReference);
}
