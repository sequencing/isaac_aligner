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
 ** \file StandardSampleSheetCsvGrammar.hpp
 **
 ** regular SampleSheet.csv grammar definition
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_STANDARD_SAMPLE_SHEET_CSV_GRAMMAR_HH
#define iSAAC_DEMULTIPLEXING_STANDARD_SAMPLE_SHEET_CSV_GRAMMAR_HH

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

namespace isaac
{
namespace demultiplexing
{


namespace bs=boost::spirit;

template <typename Iterator>
struct StandardSampleSheetCsvGrammar :
    common::CsvGrammar<Iterator>,
    bs::qi::grammar<Iterator, flowcell::BarcodeMetadataList()>
{
    StandardSampleSheetCsvGrammar(const flowcell::SequencingAdapterMetadataList &defaultAdapters) :
        StandardSampleSheetCsvGrammar::base_type(start_)
    {

        using bs::ascii::char_;
        using bs::ascii::string;
        using bs::qi::_1;
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
        using bs::uint_;

        typedef common::CsvGrammar<Iterator> Csv;

        name_ = Csv::field_;
        header_ = name_ >> *(omit[Csv::comma_] >> name_);
        header_.name("header_");

        flowcell_id_ = Csv::field_;
        flowcell_id_.name("flowcell_id_");

        lane_number_ = uint_;
        lane_number_.name("lane_number_");

        sample_name_  = Csv::field_[_val = bind(&checkIllegalCharacters, "SampleID", _1)];
        sample_name_.name("sample_name_");

        reference_ = Csv::field_;
        reference_.name("reference_");

        unknown_ = string("Undetermined")[_val = "unknown"] | //"Undetermined" is CASAVA way of specifying the sequence that does not match any known barcode
                        string("unknown")[_val = "unknown"];       //"unknown" is the recommended word to use in the sample sheet
        unknown_.name("unknown_");

        barcode_char_ = char_('-')[_val = '-'] |
                        char_('A')[_val = 'A'] | char_('a')[_val = 'A'] |
                        char_('C')[_val = 'C'] | char_('c')[_val = 'C'] |
                        char_('G')[_val = 'G'] | char_('g')[_val = 'G'] |
                        char_('T')[_val = 'T'] | char_('t')[_val = 'T'] |
                        char_('N')[_val = 'N'] | char_('n')[_val = 'N'];
        barcode_char_.name("barcode_char_");

        barcode_sequence_ = *barcode_char_;
        barcode_sequence_.name("barcode_sequence_");

//        barcode_ = unknown_ | barcode_sequence_;
        description_ = Csv::field_;
        description_.name("description_");

        control_ = char_('y')[_val = true] | char_('Y')[_val = true] | char_('n')[_val = false] | char_('N')[_val = false];
        control_.name("control_");

        recipe_ = Csv::field_;
        recipe_.name("recipe_");

        operator_ = Csv::field_;
        operator_.name("operator_");

        project_ = Csv::field_[_val = bind(&checkIllegalCharacters, "SampleProject", _1)];
        project_.name("project_");

        // unlike the real csv, sample sheets have a special treatment of commented lines
        comment_line_ = char_('#') >> *(char_(0, 0x09) | char_(0x0b, 0x0c) | char_(0x0e, 126)) >> Csv::crlf_;

        barcode_metadata_ =
            flowcell_id_[bind(&flowcell::BarcodeMetadata::setAdapters, _val, defaultAdapters)]
                         [bind(&flowcell::BarcodeMetadata::setFlowcellId, _val, _1)] >> omit[Csv::comma_] >
            lane_number_[bind(&flowcell::BarcodeMetadata::setLane, _val, _1)] > omit[Csv::comma_] >
            sample_name_[bind(&flowcell::BarcodeMetadata::setSampleName, _val, _1)] > omit[Csv::comma_] >
            reference_[bind(&flowcell::BarcodeMetadata::setReference, _val, _1)] > omit[Csv::comma_] >
            (unknown_[bind(&flowcell::BarcodeMetadata::setUnknown, _val)] |
             barcode_sequence_[bind(&flowcell::BarcodeMetadata::setSequence, _val, _1)]) > omit[Csv::comma_] >
            description_[bind(&flowcell::BarcodeMetadata::setDescription, _val, _1)] > omit[Csv::comma_] >
            control_[bind(&flowcell::BarcodeMetadata::setControl, _val, _1)] > omit[Csv::comma_] >
            recipe_[bind(&flowcell::BarcodeMetadata::setRecipe, _val, _1)] > omit[Csv::comma_] >
            operator_[bind(&flowcell::BarcodeMetadata::setOperator, _val, _1)] > omit[Csv::comma_] >
            project_[bind(&flowcell::BarcodeMetadata::setProject, _val, _1)] ;
        barcode_metadata_.name("barcode_metadata_");

        start_ = header_ >> *(+Csv::crlf_ >> *comment_line_ >> -barcode_metadata_) >> *Csv::crlf_;

    }

    bs::qi::rule<Iterator, bs::unused_type()> header_;
    bs::qi::rule<Iterator, bs::unused_type()> comment_line_;
    bs::qi::rule<Iterator, bs::unused_type()> name_;

    bs::qi::rule<Iterator, std::string()> unknown_;
    bs::qi::rule<Iterator, char()> barcode_char_;
    bs::qi::rule<Iterator, std::string()> barcode_sequence_;

    bs::qi::rule<Iterator, std::string()> project_;
    bs::qi::rule<Iterator, std::string()> operator_;
    bs::qi::rule<Iterator, std::string()> recipe_;
    bs::qi::rule<Iterator, bool()> control_;
    bs::qi::rule<Iterator, std::string()> description_;
//    bs::qi::rule<Iterator, std::string()> barcode_;
    bs::qi::rule<Iterator, std::string()> reference_;
    bs::qi::rule<Iterator, std::string()> sample_name_;
    bs::qi::rule<Iterator, unsigned()> lane_number_;
    bs::qi::rule<Iterator, std::string()> flowcell_id_;
    bs::qi::rule<Iterator, flowcell::BarcodeMetadata()> barcode_metadata_;
    bs::qi::rule<Iterator, flowcell::BarcodeMetadataList()> start_;
};

} // namespace demultiplexing
} // namespace isaac

#endif //iSAAC_DEMULTIPLEXING_STANDARD_SAMPLE_SHEET_CSV_GRAMMAR_HH
