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
 ** \file MiSeqDualBarcodeSampleSheetCsvGrammar.hpp
 **
 ** MiSeq SampleSheet.csv grammar definition for dual barcode case
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_MISEQ_DUAL_BARCODE_SAMPLE_SHEET_CSV_GRAMMAR_HH
#define iSAAC_DEMULTIPLEXING_MISEQ_DUAL_BARCODE_SAMPLE_SHEET_CSV_GRAMMAR_HH

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_attr_cast.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix1_statements.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/spirit/home/phoenix/statement/if.hpp>
#include <boost/spirit/home/phoenix/statement/throw.hpp>
#include <boost/spirit/home/phoenix/object/static_cast.hpp>
#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/spirit/home/phoenix/bind.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/support/unused.hpp>

#include "common/CsvGrammar.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "SampleSheetConstraints.hh"
#include "MiSeqSampleSheetCsvGrammar.hpp"

namespace isaac
{
namespace demultiplexing
{

namespace bs=boost::spirit;

template <typename Iterator>
struct MiSeqDualBarcodeSampleSheetCsvGrammar :
    MiSeqSampleSheetCsvGrammar<Iterator>,
    bs::qi::grammar<Iterator, flowcell::BarcodeMetadataList()>
{
        MiSeqDualBarcodeSampleSheetCsvGrammar(const flowcell::SequencingAdapterMetadataList &defaultAdapters) :
            MiSeqDualBarcodeSampleSheetCsvGrammar::base_type(start_)
    {

        using bs::ascii::char_;
        using bs::ascii::string;
        using bs::qi::_1;
        using bs::qi::_5;
        using bs::qi::omit;
        using bs::qi::_val;
        using bs::qi::eps;
//        using bs::uint_;
//        using boost::phoenix::at;
//        using boost::phoenix::back;
        using boost::phoenix::begin;
        using boost::phoenix::bind;
        using boost::phoenix::construct;
        using boost::phoenix::end;
//        using boost::phoenix::if_;
//        using boost::phoenix::insert;
        using boost::phoenix::push_back;
        using boost::phoenix::ref;
        using boost::phoenix::size;

        typedef common::CsvGrammar<Iterator> Csv;
        typedef MiSeqSampleSheetCsvGrammar<Iterator> MiSeqCsv;

        project_ = Csv::field_[_val = bind(&checkIllegalCharacters, "Sample_Project", _1)];

        column_header_ = string(
            "Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Sample_Project,index,I7_Index_ID,index2,I5_Index_ID,Description,GenomeFolder") >>
            -string(",Manifest") >> *Csv::comma_;
        column_header_.name("column_header_");

        barcode_ = MiSeqCsv::barcode_sequence_[_val = _1] > omit[Csv::comma_ > Csv::field_ > Csv::comma_] >
                   MiSeqCsv::barcode_sequence_[_val += "-"][_val += _1] > omit[Csv::comma_ > Csv::field_];
        barcode_.name("barcode_");

        barcode_metadata_ =
            MiSeqCsv::sample_id_[bind(&flowcell::BarcodeMetadata::setAdapters, _val, defaultAdapters)]
                                 [bind(&flowcell::BarcodeMetadata::setLane, _val, 1)]
                                  [bind(&flowcell::BarcodeMetadata::setOperator, _val, ref(MiSeqCsv::operatorName_))]
                                   [bind(&flowcell::BarcodeMetadata::setSampleName, _val, _1)] >> omit[Csv::comma_] >
            omit[Csv::field_] > omit[Csv::comma_] > //Sample_Name
            omit[Csv::field_] > omit[Csv::comma_] > //Sample_Plate
            omit[Csv::field_] > omit[Csv::comma_] > //Sample_Well
            project_[bind(&flowcell::BarcodeMetadata::setProject, _val, _1)] > omit[Csv::comma_] > //Sample_Project
            barcode_[bind(&flowcell::BarcodeMetadata::setSequence, _val, _1)] > omit[Csv::comma_] >
            MiSeqCsv::description_[bind(&flowcell::BarcodeMetadata::setDescription, _val, _1)] > omit[Csv::comma_] > //Description
            MiSeqCsv::reference_[bind(&flowcell::BarcodeMetadata::setReference, _val, _1)] >  //GenomeFolder
            -(omit[Csv::comma_] > omit[Csv::field_]); // optional Manifest
        barcode_metadata_.name("barcode_metadata_");

        table_ = column_header_ >> *(+Csv::crlf_ >> *MiSeqCsv::comment_line_ >> barcode_metadata_) >> *Csv::crlf_;

        data_section_ = omit[MiSeqCsv::data_section_heading_ > Csv::crlf_] >> table_;

        start_ = omit[MiSeqCsv::header_section_ > MiSeqCsv::reads_section_ >
                      -MiSeqCsv::manifests_section_ > MiSeqCsv::settings_section_] >> data_section_;
    }

    // sample sheet declarations
    bs::qi::rule<Iterator, std::string()> project_;

    bs::qi::rule<Iterator, bs::unused_type()> column_header_;
    bs::qi::rule<Iterator, std::string()> barcode_;

    bs::qi::rule<Iterator, flowcell::BarcodeMetadata()> barcode_metadata_;
    bs::qi::rule<Iterator, flowcell::BarcodeMetadataList()> data_section_;

    bs::qi::rule<Iterator, flowcell::BarcodeMetadataList()> table_;
    bs::qi::rule<Iterator, flowcell::BarcodeMetadataList()> start_;
};

} // namespace demultiplexing
} // namespace isaac

#endif //iSAAC_DEMULTIPLEXING_MISEQ_DUAL_BARCODE_SAMPLE_SHEET_CSV_GRAMMAR_HH

