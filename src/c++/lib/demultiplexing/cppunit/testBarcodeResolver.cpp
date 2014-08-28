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
 **
 ** \file testSampleSheetCsvGrammar.cpp
 **
 ** tests boost::spirit grammar for use bases mask parsing
 **
 ** \author Roman Petrovski
 **/

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testBarcodeResolver.hh"

#include "demultiplexing/BarcodeResolver.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestBarcodeResolver, registryName("BarcodeResolver"));

void TestBarcodeResolver::setUp()
{
}

void TestBarcodeResolver::tearDown()
{
}


using namespace isaac::demultiplexing;
void TestBarcodeResolver::testGet1MismatchKmer()
{
//    using isaac::demultiplexing::BarcodeResolver;
    CPPUNIT_ASSERT_EQUAL(5U * 4U, BarcodeResolver::getMismatchKmersCount(4, 1));
    CPPUNIT_ASSERT_EQUAL((5U * 4U) * (5U * 4U), BarcodeResolver::getMismatchKmersCount(4, 2));

    Kmer original = 0;
    CPPUNIT_ASSERT_EQUAL(Kmer(0), BarcodeResolver::get1MismatchKmer(original, 4, 0, 0).first);
    CPPUNIT_ASSERT_EQUAL(0U, BarcodeResolver::get1MismatchKmer(original, 4, 0, 0).second);

    CPPUNIT_ASSERT_EQUAL(Kmer(0x4), BarcodeResolver::get1MismatchKmer(original, 4, 0, 4).first);
    CPPUNIT_ASSERT_EQUAL(1U, BarcodeResolver::get1MismatchKmer(original, 4, 0, 4).second);

    CPPUNIT_ASSERT_EQUAL(Kmer(0x0), BarcodeResolver::get1MismatchKmer(original, 4, 0, 5).first);
    CPPUNIT_ASSERT_EQUAL(0U, BarcodeResolver::get1MismatchKmer(original, 4, 0, 5).second);

    CPPUNIT_ASSERT_EQUAL(Kmer(0x1 << BITS_PER_BASE), BarcodeResolver::get1MismatchKmer(original, 4, 0, 6).first);
    CPPUNIT_ASSERT_EQUAL(1U, BarcodeResolver::get1MismatchKmer(original, 4, 0, 6).second);
}

void TestBarcodeResolver::testOneComponent()
{
    std::vector<Barcode> mismatchBarcodes;
    isaac::flowcell::BarcodeMetadata barcodeMetadata;
    barcodeMetadata.setSequence("AAAA");
    barcodeMetadata.setIndex(0);
    barcodeMetadata.setComponentMismatches(std::vector<unsigned>(1, 1));
    BarcodeResolver::generateBarcodeMismatches(barcodeMetadata, mismatchBarcodes);
    CPPUNIT_ASSERT_EQUAL(Kmer(0x1 << (3 * BITS_PER_BASE)), mismatchBarcodes.at(16).getSequence());
}

std::vector<Barcode> variationsOfAdashAA = boost::assign::list_of
(Barcode(0x0000000000000000/*(AAAAAAAAAAAAAAAAAAAAA)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000001/*(AAAAAAAAAAAAAAAAAAAAC)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000002/*(AAAAAAAAAAAAAAAAAAAAG)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000003/*(AAAAAAAAAAAAAAAAAAAAT)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000004/*(AAAAAAAAAAAAAAAAAAAAN)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000008/*(AAAAAAAAAAAAAAAAAAACA)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000009/*(AAAAAAAAAAAAAAAAAAACC)*/,BarcodeId(0,0,0,0)))
(Barcode(0x000000000000000a/*(AAAAAAAAAAAAAAAAAAACG)*/,BarcodeId(0,0,0,0)))
(Barcode(0x000000000000000b/*(AAAAAAAAAAAAAAAAAAACT)*/,BarcodeId(0,0,0,0)))
(Barcode(0x000000000000000c/*(AAAAAAAAAAAAAAAAAAACN)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000010/*(AAAAAAAAAAAAAAAAAAAGA)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000011/*(AAAAAAAAAAAAAAAAAAAGC)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000012/*(AAAAAAAAAAAAAAAAAAAGG)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000013/*(AAAAAAAAAAAAAAAAAAAGT)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000014/*(AAAAAAAAAAAAAAAAAAAGN)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000018/*(AAAAAAAAAAAAAAAAAAATA)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000019/*(AAAAAAAAAAAAAAAAAAATC)*/,BarcodeId(0,0,0,0)))
(Barcode(0x000000000000001a/*(AAAAAAAAAAAAAAAAAAATG)*/,BarcodeId(0,0,0,0)))
(Barcode(0x000000000000001b/*(AAAAAAAAAAAAAAAAAAATT)*/,BarcodeId(0,0,0,0)))
(Barcode(0x000000000000001c/*(AAAAAAAAAAAAAAAAAAATN)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000020/*(AAAAAAAAAAAAAAAAAAANA)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000021/*(AAAAAAAAAAAAAAAAAAANC)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000022/*(AAAAAAAAAAAAAAAAAAANG)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000023/*(AAAAAAAAAAAAAAAAAAANT)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000024/*(AAAAAAAAAAAAAAAAAAANN)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000000/*(AAAAAAAAAAAAAAAAAAAAA)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000001/*(AAAAAAAAAAAAAAAAAAAAC)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000002/*(AAAAAAAAAAAAAAAAAAAAG)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000003/*(AAAAAAAAAAAAAAAAAAAAT)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000004/*(AAAAAAAAAAAAAAAAAAAAN)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000040/*(AAAAAAAAAAAAAAAAAACAA)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000041/*(AAAAAAAAAAAAAAAAAACAC)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000042/*(AAAAAAAAAAAAAAAAAACAG)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000043/*(AAAAAAAAAAAAAAAAAACAT)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000044/*(AAAAAAAAAAAAAAAAAACAN)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000080/*(AAAAAAAAAAAAAAAAAAGAA)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000081/*(AAAAAAAAAAAAAAAAAAGAC)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000082/*(AAAAAAAAAAAAAAAAAAGAG)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000083/*(AAAAAAAAAAAAAAAAAAGAT)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000084/*(AAAAAAAAAAAAAAAAAAGAN)*/,BarcodeId(0,0,0,0)))
(Barcode(0x00000000000000c0/*(AAAAAAAAAAAAAAAAAATAA)*/,BarcodeId(0,0,0,0)))
(Barcode(0x00000000000000c1/*(AAAAAAAAAAAAAAAAAATAC)*/,BarcodeId(0,0,0,0)))
(Barcode(0x00000000000000c2/*(AAAAAAAAAAAAAAAAAATAG)*/,BarcodeId(0,0,0,0)))
(Barcode(0x00000000000000c3/*(AAAAAAAAAAAAAAAAAATAT)*/,BarcodeId(0,0,0,0)))
(Barcode(0x00000000000000c4/*(AAAAAAAAAAAAAAAAAATAN)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000100/*(AAAAAAAAAAAAAAAAAANAA)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000101/*(AAAAAAAAAAAAAAAAAANAC)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000102/*(AAAAAAAAAAAAAAAAAANAG)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000103/*(AAAAAAAAAAAAAAAAAANAT)*/,BarcodeId(0,0,0,0)))
(Barcode(0x0000000000000104/*(AAAAAAAAAAAAAAAAAANAN)*/,BarcodeId(0,0,0,0)))
;

void TestBarcodeResolver::testTwoComponents()
{
    std::vector<Barcode> mismatchBarcodes;
    isaac::flowcell::BarcodeMetadata barcodeMetadata;
    barcodeMetadata.setSequence("AA-A");
    barcodeMetadata.setIndex(0);
    barcodeMetadata.setComponentMismatches(std::vector<unsigned>(2, 1));
    BarcodeResolver::generateBarcodeMismatches(barcodeMetadata, mismatchBarcodes);

    BOOST_FOREACH(const Barcode &test, variationsOfAdashAA)
    {
        CPPUNIT_ASSERT_EQUAL(test.getSequence(),
                             mismatchBarcodes.at(&test - &variationsOfAdashAA.front()).getSequence());
    }
}

void TestBarcodeResolver::testMismatchCollision()
{
    isaac::flowcell::BarcodeMetadataList barcodeMetadataList(3);
    std::vector<unsigned> compMism(2, 1);
    barcodeMetadataList.at(0).setUnknown();
    barcodeMetadataList.at(0).setIndex(0);
    barcodeMetadataList.at(0).setComponentMismatches(compMism);
    barcodeMetadataList.at(1).setSequence("G-AA");
    barcodeMetadataList.at(1).setIndex(1);
    barcodeMetadataList.at(1).setComponentMismatches(compMism);
    barcodeMetadataList.at(2).setSequence("T-CC");
    barcodeMetadataList.at(2).setIndex(2);
    barcodeMetadataList.at(2).setComponentMismatches(compMism);

    CPPUNIT_ASSERT_THROW(BarcodeResolver::generateMismatches(barcodeMetadataList, barcodeMetadataList),
                         isaac::common::InvalidOptionException);

}

