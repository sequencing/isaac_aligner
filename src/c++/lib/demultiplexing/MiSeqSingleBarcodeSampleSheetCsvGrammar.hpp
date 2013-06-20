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
 ** \file MiSeqSingleBarcodeSampleSheetCsvGrammar.hpp
 **
 ** MiSeq SampleSheet.csv grammar definition for dual barcode case
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_MISEQ_SINGLE_BARCODE_SAMPLE_SHEET_CSV_GRAMMAR_HH
#define iSAAC_DEMULTIPLEXING_MISEQ_SINGLE_BARCODE_SAMPLE_SHEET_CSV_GRAMMAR_HH

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
struct MiSeqSingleBarcodeSampleSheetCsvGrammar :
    MiSeqSampleSheetCsvGrammar<Iterator>,
    bs::qi::grammar<Iterator, flowcell::BarcodeMetadataList()>
{
        MiSeqSingleBarcodeSampleSheetCsvGrammar(const flowcell::SequencingAdapterMetadataList &defaultAdapters) :
            MiSeqSingleBarcodeSampleSheetCsvGrammar::base_type(start_)
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

        // column header determines the type of the sample sheet. All subsequent rules use expectation parsers wherever possible
        column_header_ = string("Sample_ID,Sample_Name,GenomeFolder,Index") >> -string(",Manifest") >> *Csv::comma_;

        barcode_metadata_ =
            MiSeqCsv::sample_id_[bind(&flowcell::BarcodeMetadata::setAdapters, _val, defaultAdapters)]
                                 [bind(&flowcell::BarcodeMetadata::setLane, _val, 1)]
                                  [bind(&flowcell::BarcodeMetadata::setOperator, _val, ref(MiSeqCsv::operatorName_))]
                                   [bind(&flowcell::BarcodeMetadata::setProject, _val, ref(MiSeqCsv::projectName_))]
                                    [bind(&flowcell::BarcodeMetadata::setSampleName, _val, _1)] >> omit[Csv::comma_] > //Sample_ID
            omit[Csv::field_] > omit[Csv::comma_] > //Sample_Name
            MiSeqCsv::reference_[bind(&flowcell::BarcodeMetadata::setReference, _val, _1)]  > omit[Csv::comma_] > //GenomeFolder
            MiSeqCsv::barcode_sequence_[bind(&flowcell::BarcodeMetadata::setSequence, _val, _1)] > //Index
            -(omit[Csv::comma_] > omit[Csv::field_]); // optional Manifest

        table_ = column_header_ >> *(+Csv::crlf_ >> *MiSeqCsv::comment_line_ >> barcode_metadata_) >> *Csv::crlf_;

        data_section_ = omit[MiSeqCsv::data_section_heading_ > Csv::crlf_] >> table_;

        start_ = omit[MiSeqCsv::header_section_ > MiSeqCsv::reads_section_ >
                      -MiSeqCsv::manifests_section_ > MiSeqCsv::settings_section_] >> data_section_;
    }

    // sample sheet declarations
    bs::qi::rule<Iterator, bs::unused_type()> column_header_;

    bs::qi::rule<Iterator, flowcell::BarcodeMetadata()> barcode_metadata_;
    bs::qi::rule<Iterator, flowcell::BarcodeMetadataList()> data_section_;

    bs::qi::rule<Iterator, flowcell::BarcodeMetadataList()> table_;
    bs::qi::rule<Iterator, flowcell::BarcodeMetadataList()> start_;
};

} // namespace demultiplexing
} // namespace isaac

#endif //iSAAC_DEMULTIPLEXING_MISEQ_SINGLE_BARCODE_SAMPLE_SHEET_CSV_GRAMMAR_HH

