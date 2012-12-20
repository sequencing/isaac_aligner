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
#include <vector>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testSortedReferenceXml.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestSortedReferenceXml, registryName("SortedReferenceXml"));

TestSortedReferenceXml::TestSortedReferenceXml()
    : xmlString(
"<?xml version=\"1.0\"?>\n"
"<SortedReference>\n"
"  <FormatVersion>2</FormatVersion>\n"
"  <Contigs>\n"
"    <Contig Position=\"0\">\n"
"      <Index>0</Index>\n"
"      <Name>10</Name>\n"
"      <Sequence>\n"
"        <File>/illumina/scratch/smallprojects/iSAAC/Genomes/HumanNCBI36Fasta_all.fa</File>\n"
"        <Offset>54</Offset>\n"
"        <Size>137630983</Size>\n"
"      </Sequence>\n"
"      <TotalBases>135374737</TotalBases>\n"
"      <AcgtBases>131624728</AcgtBases>\n"
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

void TestSortedReferenceXml::testGetDefaultMaskWidth()
{
    std::istringstream is(xmlString);
    isaac::reference::SortedReferenceXml sortedReferenceXml;
    is >> sortedReferenceXml;
    CPPUNIT_ASSERT_EQUAL(1U, sortedReferenceXml.getDefaultMaskWidth());
}

void TestSortedReferenceXml::testGetMaxPrefixRangeCount()
{
    std::istringstream is(xmlString);
    isaac::reference::SortedReferenceXml sortedReferenceXml;
    is >> sortedReferenceXml;
    CPPUNIT_ASSERT_EQUAL(927297U, sortedReferenceXml.getMaxPrefixRangeCount());
}

void TestSortedReferenceXml::testGetMaskFileList()
{
    std::istringstream is(xmlString);
    isaac::reference::SortedReferenceXml sortedReferenceXml;
    is >> sortedReferenceXml;
    std::vector<isaac::reference::SortedReferenceXml::MaskFile> list =
        sortedReferenceXml.getMaskFileList("ABCD");
    CPPUNIT_ASSERT_EQUAL(3U, unsigned(list.size()));
    CPPUNIT_ASSERT_EQUAL(boost::filesystem::path("/illumina/scratch/smallprojects/iSAAC/Genomes/HumanNCBI36-neighbors/HumanNCBI36Fasta_all.fa-32mer-6bit-ABCD-02.dat"),
                         list.back().path);
    CPPUNIT_ASSERT_EQUAL(sortedReferenceXml.getDefaultMaskWidth(), list.back().maskWidth);
}
