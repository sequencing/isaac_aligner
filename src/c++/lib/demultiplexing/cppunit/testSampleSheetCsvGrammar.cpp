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
#include "testSampleSheetCsvGrammar.hh"

#include "SampleSheetCsvGrammar.hpp"
#include "flowcell/SequencingAdapterMetadata.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestSampleSheetCsvGrammar, registryName("SampleSheetCsvGrammar"));

void TestSampleSheetCsvGrammar::setUp()
{
}

void TestSampleSheetCsvGrammar::tearDown()
{
}


using isaac::demultiplexing::SampleSheetCsvGrammar;
typedef SampleSheetCsvGrammar<std::string::const_iterator> Parser;

void TestSampleSheetCsvGrammar::testStandard()
{
    std::string test(
        // header line (first line) must be always discarded
        "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\r"
        // commented lines must be discarded
        "#A805CKABXX,1,AR008,human,ACTTGA,Cypress,Y,\"101,7,101\",CB,Demo\r\n"
        // regular line ending with just lf
        "A805CKABXX,1,AR005,human,ACTTGA,Cypress,Y,101+7,CB,Demo\n"
        // regular line ending with crlf and containing quoted fields with quotes and comma inside
        "A805CKABXX,2,AR008,human,ACT-TGA,\"\"\"Cypress\"\"\",Y,\"101,7,101\",CB,Demo\r\n"
        // empty line
        "\r\n"
        // another empty line
        "\r"
        // yet another empty line
        "\n"
        // another regular line ending with cr only
        "A805CKABXX,3,AR008,human,ACTTGA,Cypress,Y,101+7,CB,Demo\r"
        "A805CKABXX,4,,human,unknown,Cypress,Y,101+7,CB,Demo\r"
        "A805CKABXX,5,,human,Undetermined,Cypress,Y,101+7,CB,Demo\r"
        "#A805CKABXX,6,,human,Undetermined,Cypress,Y,101+7,CB,Demo\r"
        "#A805CKABXX,7,,human,Undetermined,Cypress,Y,101+7,CB,Demo\r"
        );
    std::string::const_iterator parseIt(test.begin());
    std::string::const_iterator parseEnd(test.end());

    const isaac::flowcell::SequencingAdapterMetadata testAdapterMetadata(
        "CTGTCTCTTATACACATCT",
        false,
        strlen("CTGTCTCTTATACACATCT"));

    static const isaac::flowcell::SequencingAdapterMetadataList adapters = boost::assign::list_of(testAdapterMetadata);
    SampleSheetCsvGrammar<std::string::const_iterator> parser(adapters);
    isaac::flowcell::BarcodeMetadataList result;

    try
    {
        if (!boost::spirit::qi::parse(parseIt, parseEnd, parser, result) || parseEnd != parseIt)
        {
            CPPUNIT_FAIL("Could not parse the sample sheet csv stream text:\n" + std::string(parseIt, parseEnd));
        }
    }
    catch(boost::spirit::qi::expectation_failure<std::string::const_iterator> const &e)
    {
        const std::string::const_iterator bufferBegin(test.begin());

        const std::string message = (boost::format("Could not parse the sample sheet csv. "
            "Expected:"
            "\n==================================================================\n%s||"
            "\n==================================================================\n"
            "Got:"
            "\n==================================================================\n%s||"
            "\n==================================================================\n"
            " at offset: %d") % e.what_ % std::string(e.first, e.last) %
            std::distance(bufferBegin, e.first)).str();
        CPPUNIT_FAIL(message);
    }

    CPPUNIT_ASSERT_EQUAL(size_t(5), result.size());

    CPPUNIT_ASSERT_EQUAL(std::string("A805CKABXX"), result.at(0).getFlowcellId());
    CPPUNIT_ASSERT_EQUAL(1U, result.at(0).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("AR005"), result.at(0).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("ACTTGA"), result.at(0).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("Cypress"), result.at(0).getDescription());
    CPPUNIT_ASSERT_EQUAL(std::string("101+7"), result.at(0).getRecipe());
    CPPUNIT_ASSERT_EQUAL(std::string("CB"), result.at(0).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Demo"), result.at(0).getProject());
    CPPUNIT_ASSERT_EQUAL(size_t(1), result.at(0).getAdapters().size());
    CPPUNIT_ASSERT_EQUAL(testAdapterMetadata, result.at(0).getAdapters().at(0));
    CPPUNIT_ASSERT_EQUAL(false, result.at(0).isUnknown());

    CPPUNIT_ASSERT_EQUAL(std::string("A805CKABXX"), result.at(1).getFlowcellId());
    CPPUNIT_ASSERT_EQUAL(2U, result.at(1).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("AR008"), result.at(1).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("ACT-TGA"), result.at(1).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("\"Cypress\""), result.at(1).getDescription());
    CPPUNIT_ASSERT_EQUAL(std::string("101,7,101"), result.at(1).getRecipe());
    CPPUNIT_ASSERT_EQUAL(std::string("CB"), result.at(1).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Demo"), result.at(1).getProject());
    CPPUNIT_ASSERT_EQUAL(size_t(1), result.at(1).getAdapters().size());
    CPPUNIT_ASSERT_EQUAL(testAdapterMetadata, result.at(1).getAdapters().at(0));
    CPPUNIT_ASSERT_EQUAL(false, result.at(1).isUnknown());

    CPPUNIT_ASSERT_EQUAL(std::string("A805CKABXX"), result.at(2).getFlowcellId());
    CPPUNIT_ASSERT_EQUAL(3U, result.at(2).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("AR008"), result.at(2).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("ACTTGA"), result.at(2).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("Cypress"), result.at(2).getDescription());
    CPPUNIT_ASSERT_EQUAL(std::string("101+7"), result.at(2).getRecipe());
    CPPUNIT_ASSERT_EQUAL(std::string("CB"), result.at(2).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Demo"), result.at(2).getProject());
    CPPUNIT_ASSERT_EQUAL(size_t(1), result.at(2).getAdapters().size());
    CPPUNIT_ASSERT_EQUAL(testAdapterMetadata, result.at(2).getAdapters().at(0));
    CPPUNIT_ASSERT_EQUAL(false, result.at(2).isUnknown());

    CPPUNIT_ASSERT_EQUAL(std::string("A805CKABXX"), result.at(3).getFlowcellId());
    CPPUNIT_ASSERT_EQUAL(4U, result.at(3).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("unknown"), result.at(3).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string(""), result.at(3).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("Cypress"), result.at(3).getDescription());
    CPPUNIT_ASSERT_EQUAL(std::string("101+7"), result.at(3).getRecipe());
    CPPUNIT_ASSERT_EQUAL(std::string("CB"), result.at(3).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Demo"), result.at(3).getProject());
    CPPUNIT_ASSERT_EQUAL(size_t(1), result.at(3).getAdapters().size());
    CPPUNIT_ASSERT_EQUAL(testAdapterMetadata, result.at(3).getAdapters().at(0));
    CPPUNIT_ASSERT_EQUAL(true, result.at(3).isUnknown());

    CPPUNIT_ASSERT_EQUAL(std::string("A805CKABXX"), result.at(4).getFlowcellId());
    CPPUNIT_ASSERT_EQUAL(5U, result.at(4).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("unknown"), result.at(4).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string(""), result.at(4).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("Cypress"), result.at(4).getDescription());
    CPPUNIT_ASSERT_EQUAL(std::string("101+7"), result.at(4).getRecipe());
    CPPUNIT_ASSERT_EQUAL(std::string("CB"), result.at(4).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Demo"), result.at(4).getProject());
    CPPUNIT_ASSERT_EQUAL(size_t(1), result.at(4).getAdapters().size());
    CPPUNIT_ASSERT_EQUAL(testAdapterMetadata, result.at(4).getAdapters().at(0));
    CPPUNIT_ASSERT_EQUAL(true, result.at(4).isUnknown());
}

void TestSampleSheetCsvGrammar::testDualBarcodeMiSeq()
{
    std::string test(
        "[Header],,,,,,,,,,\r"
"IEMFileVersion,3,,,,,,,,,\n"
"Investigator Name,Isabelle,,,,,,,,,\n"
"Project Name,Zebra_validation,,,,,,,,,\n"
"Experiment Name,48plex,,,,,,,,,\n"
"Date,20/02/2012,,,,,,,,,\n"
"Workflow,Resequencing,,,,,,,,,\n"
"Assay,TruSeq,,,,,,,,,\n"
"Description,G7_H1_12pM,,,,,,,,,\n"
"Chemistry,Amplicon,,,,,,,,,\n"
"[Reads],,,,,,,,,,\n"
"151,,,,,,,,,,\n"
"151,,,,,,,,,,\n"
"[Settings],,,,,,,,,,\n"
"[Data],,,,,,,,,,\n"
"Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Sample_Project,index,I7_Index_ID,index2,I5_Index_ID,Description,GenomeFolder\n"
"A1_Bcereus,Human,A1_Bcereus,H12,Zebra_validation,ATTACTCG,D701,TATAGCCT,D501,Tada,\\\\ch-isilon\\isilon\\Mondas_software\\Genomes\\B_Cereus_ATCC10987_ELAND\n"
"C2_Bcereus,Human,C2_Bcereus,H11,Zebra_validation,TCCGGAGA,D702,CCTATCCT,D503,none,\\\\ch-isilon\\isilon\\Mondas_software\\Genomes\\B_Cereus_ATCC10987_ELAND\n"
        );
    std::string::const_iterator parseIt(test.begin());
    std::string::const_iterator parseEnd(test.end());

    const isaac::flowcell::SequencingAdapterMetadataList noAdapters;
    SampleSheetCsvGrammar<std::string::const_iterator> parser(noAdapters);
    isaac::flowcell::BarcodeMetadataList result;
    try
    {
        if (!boost::spirit::qi::parse(parseIt, parseEnd, parser, result) || parseEnd != parseIt)
        {
            CPPUNIT_FAIL("Could not parse the sample sheet csv stream text:\n" + std::string(parseIt, parseEnd));
        }
    }
    catch(boost::spirit::qi::expectation_failure<std::string::const_iterator> const &e)
    {
        const std::string::const_iterator bufferBegin(test.begin());

        const std::string message = (boost::format("Could not parse the sample sheet csv. "
            "Expected:"
            "\n==================================================================\n%s||"
            "\n==================================================================\n"
            "Got:"
            "\n==================================================================\n%s||"
            "\n==================================================================\n"
            " at offset: %d") % e.what_ % std::string(e.first, e.last) %
            std::distance(bufferBegin, e.first)).str();
        CPPUNIT_FAIL(message);
    }
    CPPUNIT_ASSERT_EQUAL(size_t(2), result.size());

    CPPUNIT_ASSERT_EQUAL(1U, result.at(0).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("A1_Bcereus"), result.at(0).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("ATTACTCG-TATAGCCT"), result.at(0).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("Tada"), result.at(0).getDescription());
    CPPUNIT_ASSERT_EQUAL(std::string("Isabelle"), result.at(0).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Zebra_validation"), result.at(0).getProject());
    CPPUNIT_ASSERT_EQUAL(false, result.at(0).isUnknown());

    CPPUNIT_ASSERT_EQUAL(1U, result.at(1).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("C2_Bcereus"), result.at(1).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("TCCGGAGA-CCTATCCT"), result.at(1).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("none"), result.at(1).getDescription());
    CPPUNIT_ASSERT_EQUAL(std::string("Isabelle"), result.at(1).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Zebra_validation"), result.at(1).getProject());
    CPPUNIT_ASSERT_EQUAL(false, result.at(1).isUnknown());
}

void TestSampleSheetCsvGrammar::testAnotherDualBarcodeMiSeq()
{
    std::string test(
        "[Header],\n\r"
"IEMFileVersion,3\n\r"
"Investigator Name,Nick\n\r"
"Project Name,Tudu\n\r"
"Experiment Name,DVT_VaraiabilityHC#11_SD\n\r"
"Date,1/18/2012\n\r"
"Workflow,Resequencing\n\r"
"Assay,Nextera\n\r"
"Description,Variability\n\r"
"Chemistry,Amplicon\n\r"
"\n\r"
"[Reads],\n\r"
"151,\n\r"
"151,\n\r"
"\n\r"
"[Manifests],\n\r"
"A,Manifest,,,,,,,,,,\n\r"
"\n\r"
"[Settings],,,,,,,,,,,\n\r"
"Aligner,isaac,,,,,,,,,,\n\r"
"Adapter,CTGTCTCTTATACACATCT,,,,,,,,,,\n\r"
"\n\r"
"[Data],,,,,,,,,,,\n\r"
"Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Sample_Project,index,I7_Index_ID,index2,I5_Index_ID,Description,GenomeFolder,Manifest\n\r"
"V_1,V_1,Variability,F09,Tudu,GCTACGCT,N709,ACTGCATA,N506,DaytoDayVaraiability,Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFASTA/,A\n\r"
"V_2,V_2,Variability,F10,Tudu,CGAGGCTG,N710,ACTGCATA,N506,DaytoDayVaraiability,Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFASTA/,A\n\r"
        );
    std::string::const_iterator parseIt(test.begin());
    std::string::const_iterator parseEnd(test.end());

    const isaac::flowcell::SequencingAdapterMetadataList noAdapters;
    SampleSheetCsvGrammar<std::string::const_iterator> parser(noAdapters);
    isaac::flowcell::BarcodeMetadataList result;
    try
    {
        if (!boost::spirit::qi::parse(parseIt, parseEnd, parser, result) || parseEnd != parseIt)
        {
            CPPUNIT_FAIL("Could not parse the sample sheet csv stream text:\n" + std::string(parseIt, parseEnd));
        }
    }
    catch(boost::spirit::qi::expectation_failure<std::string::const_iterator> const &e)
    {
        const std::string::const_iterator bufferBegin(test.begin());

        const std::string message = (boost::format("Could not parse the sample sheet csv. "
            "Expected:"
            "\n==================================================================\n%s||"
            "\n==================================================================\n"
            "Got:"
            "\n==================================================================\n%s||"
            "\n==================================================================\n"
            " at offset: %d") % e.what_ % std::string(e.first, e.last) %
            std::distance(bufferBegin, e.first)).str();
        CPPUNIT_FAIL(message);
    }
    CPPUNIT_ASSERT_EQUAL(size_t(2), result.size());

    CPPUNIT_ASSERT_EQUAL(1U, result.at(0).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("V_1"), result.at(0).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("GCTACGCT-ACTGCATA"), result.at(0).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("DaytoDayVaraiability"), result.at(0).getDescription());
    CPPUNIT_ASSERT_EQUAL(std::string("Nick"), result.at(0).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Tudu"), result.at(0).getProject());
    CPPUNIT_ASSERT_EQUAL(false, result.at(0).isUnknown());

    CPPUNIT_ASSERT_EQUAL(1U, result.at(1).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("V_2"), result.at(1).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("CGAGGCTG-ACTGCATA"), result.at(1).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("DaytoDayVaraiability"), result.at(1).getDescription());
    CPPUNIT_ASSERT_EQUAL(std::string("Nick"), result.at(1).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Tudu"), result.at(1).getProject());
    CPPUNIT_ASSERT_EQUAL(false, result.at(1).isUnknown());
}


void TestSampleSheetCsvGrammar::testNonMultiplexedMiSeq()
{
    std::string test(
"[Header],,,,,,,,,,,,,,,,\r\n"
"IEMFileVersion,3,,,,,,,,,,,,,,,\r\n"
"Investigator Name,TN/ZK,,,,,,,,,,,,,,,\r\n"
"Project Name,2kb,,,,,,,,,,,,,,,\r\n"
"Experiment Name,Roberto 2kb,,,,,,,,,,,,,,,\r\n"
"Date,05/01/2012,,,,,,,,,,,,,,,\r\n"
"Workflow,Resequencing,,,,,,,,,,,,,,,\r\n"
"Assay,TruSeq DNA/RNA,,,,,,,,,,,,,,,\r\n"
"Description,,,,,,,,,,,,,,,,\r\n"
"Chemistry,Arusha_Roberto,,,,,,,,,,,,,,,\r\n"
"[Reads],,,,,,,,,,,,,,,,\r\n"
"101,,,,,,,,,,,,,,,,\r\n"
"101,,,,,,,,,,,,,,,,\r\n"
"[Settings],,,,,,,,,,,,,,,,\r\n"
"[Data],,,,,,,,,,,,,,,,\r\n"
"Sample_ID,Sample_Name,GenomeFolder\r\n"
"E-coli,CT5244,\\\\ch-isilon\\isilon\\Mondas_software\\Genomes\\E_coli_ELAND"
        );
    std::string::const_iterator parseIt(test.begin());
    std::string::const_iterator parseEnd(test.end());

    const isaac::flowcell::SequencingAdapterMetadataList noAdapters;
    SampleSheetCsvGrammar<std::string::const_iterator> parser(noAdapters);
    isaac::flowcell::BarcodeMetadataList result;
    try
    {
        if (!boost::spirit::qi::parse(parseIt, parseEnd, parser, result) || parseEnd != parseIt)
        {
            CPPUNIT_FAIL("Could not parse the sample sheet csv stream text:\n" + std::string(parseIt, parseEnd));
        }
    }
    catch(boost::spirit::qi::expectation_failure<std::string::const_iterator> const &e)
    {
        const std::string::const_iterator bufferBegin(test.begin());

        const std::string message = (boost::format("Could not parse the sample sheet csv. "
            "Expected:"
            "\n==================================================================\n%s||"
            "\n==================================================================\n"
            "Got:"
            "\n==================================================================\n%s||"
            "\n==================================================================\n"
            " at offset: %d") % e.what_ % std::string(e.first, e.last) %
            std::distance(bufferBegin, e.first)).str();
        CPPUNIT_FAIL(message);
    }
    CPPUNIT_ASSERT_EQUAL(size_t(1), result.size());

    CPPUNIT_ASSERT_EQUAL(1U, result.at(0).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("E-coli"), result.at(0).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string(""), result.at(0).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("none"), result.at(0).getName());
    CPPUNIT_ASSERT_EQUAL(std::string("TN/ZK"), result.at(0).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("2kb"), result.at(0).getProject());
    CPPUNIT_ASSERT_EQUAL(false, result.at(0).isUnknown());
    CPPUNIT_ASSERT_EQUAL(true, result.at(0).isNoIndex());
    CPPUNIT_ASSERT_EQUAL(std::string("\\\\ch-isilon\\isilon\\Mondas_software\\Genomes\\E_coli_ELAND"), result.at(0).getReference());
}


void TestSampleSheetCsvGrammar::testSingleBarcodeMiSeq()
{
    std::string test(
        "[Header],,,\r\n"
        "Investigator Name,Aurelie,,\r\n"
        "Project Name,Blah,,\r\n"
        "Experiment Name,tada,,\r\n"
        "Date,02/03/2012,,\r\n"
        "Workflow,Resequencing,,\r\n"
        "Chemistry,outch,,\r\n"
        ",,,\r\n"
        "[Reads],,,\r\n"
        "251,,,\r\n"
        "251,,,\r\n"
        ",,,\r\n"
        "[Manifests],,,\r\n"
        "A,ManifestNameHere,,\r\n"
        ",,,\r\n"
        "[Settings],,,\r\n"
        "FilterPCRDuplicates,0,,\r\n"
        ",,,\r\n"
        "[Data],,,\r\n"
        "Sample_ID,Sample_Name,GenomeFolder,Index\r\n"
        "1,BCereus,\\\\ch-isilon\\iGenomes\\Bacillus_cereus_ATCC_10987\\NCBI\\2004-02-13\\Sequence\\Chromosomes,CTTGTA\r\n"
        "2,Rhodo,\\\\ch-isilon\\iGenomes\\Rhodobacter_sphaeroides_2.4.1\\NCBI\\2005-10-07\\Sequence\\Chromosomes,CAGATC\r\n"
        "3,Human,\\\\ch-isilon\\iGenomes\\Homo_sapiens\\NCBI\\build37.2\\Sequence\\Chromosomes,ATCACG\r\n"
        "4,EColi,\\\\ch-isilon\\iGenomes\\Escherichia_coli_K_12_DH10B\\NCBI\\2008-03-17\\Sequence\\Chromosomes,TGACCA\r\n"
        "5,EColi,\\\\ch-isilon\\iGenomes\\Escherichia_coli_K_12_DH10B\\NCBI\\2008-03-17\\Sequence\\Chromosomes,GCCAAT\r\n"
        "6,PHix,\\\\ch-isilon\\iGenomes\\PhiX\\Illumina\\RTA\\Sequence\\Chromosomes,CGATGT\r\n"
        );
    std::string::const_iterator parseIt(test.begin());
    std::string::const_iterator parseEnd(test.end());

    const isaac::flowcell::SequencingAdapterMetadataList noAdapters;
    SampleSheetCsvGrammar<std::string::const_iterator> parser(noAdapters);
    isaac::flowcell::BarcodeMetadataList result;
    try
    {
        if (!boost::spirit::qi::parse(parseIt, parseEnd, parser, result) || parseEnd != parseIt)
        {
            CPPUNIT_FAIL("Could not parse the sample sheet csv stream text:\n" + std::string(parseIt, parseEnd));
        }
    }
    catch(boost::spirit::qi::expectation_failure<std::string::const_iterator> const &e)
    {
        const std::string::const_iterator bufferBegin(test.begin());

        const std::string message = (boost::format("Could not parse the sample sheet csv. "
            "Expected:"
            "\n==================================================================\n%s||"
            "\n==================================================================\n"
            "Got:"
            "\n==================================================================\n%s||"
            "\n==================================================================\n"
            " at offset: %d") % e.what_ % std::string(e.first, e.last) %
            std::distance(bufferBegin, e.first)).str();
        CPPUNIT_FAIL(message);
    }
    CPPUNIT_ASSERT_EQUAL(size_t(6), result.size());

    CPPUNIT_ASSERT_EQUAL(1U, result.at(0).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("1"), result.at(0).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("CTTGTA"), result.at(0).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("CTTGTA"), result.at(0).getName());
    CPPUNIT_ASSERT_EQUAL(std::string("Aurelie"), result.at(0).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Blah"), result.at(0).getProject());
    CPPUNIT_ASSERT_EQUAL(false, result.at(0).isUnknown());
    CPPUNIT_ASSERT_EQUAL(false, result.at(0).isNoIndex());
    CPPUNIT_ASSERT_EQUAL(std::string("\\\\ch-isilon\\iGenomes\\Bacillus_cereus_ATCC_10987\\NCBI\\2004-02-13\\Sequence\\Chromosomes"), result.at(0).getReference());

    CPPUNIT_ASSERT_EQUAL(1U, result.at(1).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("2"), result.at(1).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("CAGATC"), result.at(1).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("CAGATC"), result.at(1).getName());
    CPPUNIT_ASSERT_EQUAL(std::string("Aurelie"), result.at(1).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Blah"), result.at(1).getProject());
    CPPUNIT_ASSERT_EQUAL(false, result.at(1).isUnknown());
    CPPUNIT_ASSERT_EQUAL(false, result.at(1).isNoIndex());
    CPPUNIT_ASSERT_EQUAL(std::string("\\\\ch-isilon\\iGenomes\\Rhodobacter_sphaeroides_2.4.1\\NCBI\\2005-10-07\\Sequence\\Chromosomes"), result.at(1).getReference());

    CPPUNIT_ASSERT_EQUAL(1U, result.at(2).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("3"), result.at(2).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("ATCACG"), result.at(2).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("ATCACG"), result.at(2).getName());
    CPPUNIT_ASSERT_EQUAL(std::string("Aurelie"), result.at(2).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Blah"), result.at(2).getProject());
    CPPUNIT_ASSERT_EQUAL(false, result.at(2).isUnknown());
    CPPUNIT_ASSERT_EQUAL(false, result.at(2).isNoIndex());
    CPPUNIT_ASSERT_EQUAL(std::string("\\\\ch-isilon\\iGenomes\\Homo_sapiens\\NCBI\\build37.2\\Sequence\\Chromosomes"), result.at(2).getReference());

    CPPUNIT_ASSERT_EQUAL(1U, result.at(3).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("4"), result.at(3).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("TGACCA"), result.at(3).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("TGACCA"), result.at(3).getName());
    CPPUNIT_ASSERT_EQUAL(std::string("Aurelie"), result.at(3).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Blah"), result.at(3).getProject());
    CPPUNIT_ASSERT_EQUAL(false, result.at(3).isUnknown());
    CPPUNIT_ASSERT_EQUAL(false, result.at(3).isNoIndex());
    CPPUNIT_ASSERT_EQUAL(std::string("\\\\ch-isilon\\iGenomes\\Escherichia_coli_K_12_DH10B\\NCBI\\2008-03-17\\Sequence\\Chromosomes"), result.at(3).getReference());

    CPPUNIT_ASSERT_EQUAL(1U, result.at(4).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("5"), result.at(4).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("GCCAAT"), result.at(4).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("GCCAAT"), result.at(4).getName());
    CPPUNIT_ASSERT_EQUAL(std::string("Aurelie"), result.at(4).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Blah"), result.at(4).getProject());
    CPPUNIT_ASSERT_EQUAL(false, result.at(4).isUnknown());
    CPPUNIT_ASSERT_EQUAL(false, result.at(4).isNoIndex());
    CPPUNIT_ASSERT_EQUAL(std::string("\\\\ch-isilon\\iGenomes\\Escherichia_coli_K_12_DH10B\\NCBI\\2008-03-17\\Sequence\\Chromosomes"), result.at(4).getReference());

    CPPUNIT_ASSERT_EQUAL(1U, result.at(5).getLane());
    CPPUNIT_ASSERT_EQUAL(std::string("6"), result.at(5).getSampleName());
    CPPUNIT_ASSERT_EQUAL(std::string("CGATGT"), result.at(5).getSequence());
    CPPUNIT_ASSERT_EQUAL(std::string("CGATGT"), result.at(5).getName());
    CPPUNIT_ASSERT_EQUAL(std::string("Aurelie"), result.at(5).getOperator());
    CPPUNIT_ASSERT_EQUAL(std::string("Blah"), result.at(5).getProject());
    CPPUNIT_ASSERT_EQUAL(false, result.at(5).isUnknown());
    CPPUNIT_ASSERT_EQUAL(false, result.at(5).isNoIndex());
    CPPUNIT_ASSERT_EQUAL(std::string("\\\\ch-isilon\\iGenomes\\PhiX\\Illumina\\RTA\\Sequence\\Chromosomes"), result.at(5).getReference());
}
