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
 ** \file MiSeqSampleSheetCsvGrammar.hpp
 **
 ** MiSeq SampleSheet.csv grammar definition
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_MISEQ_SAMPLE_SHEET_CSV_GRAMMAR_HH
#define iSAAC_DEMULTIPLEXING_MISEQ_SAMPLE_SHEET_CSV_GRAMMAR_HH

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
#include "SampleSheetConstraints.hh"
#include "flowcell/BarcodeMetadata.hh"

namespace isaac
{
namespace demultiplexing
{

namespace bs=boost::spirit;

template <typename Iterator>
struct MiSeqSampleSheetCsvGrammar :
    common::CsvGrammar<Iterator>
{
    MiSeqSampleSheetCsvGrammar()
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

        header_section_heading_ = char_('[') >> string("Header") >> char_(']') >> *(omit[Csv::comma_] >> Csv::field_);

        operator_ = omit[string("Investigator Name")] >> omit[Csv::comma_] >>
            Csv::field_[ref(operatorName_) = _1] >> *(omit[Csv::comma_] >> omit[Csv::field_]);
        operator_.name("operator_");
        project_ = omit[string("Project Name")] >> omit[Csv::comma_] >>
            Csv::field_[ref(projectName_) = _1] >> *(omit[Csv::comma_] >> omit[Csv::field_]);
        project_.name("project_");

        header_section_ = header_section_heading_ >> Csv::crlf_ >>
            *((section_line_ - operator_ - project_) >> Csv::crlf_) >>
            operator_ >> Csv::crlf_ >> project_ >> Csv::crlf_ >>
            *((section_line_ - operator_ - project_) >> Csv::crlf_);
        header_section_.name("header_section_");

        reads_section_heading_ = char_('[') >> string("Reads") >> char_(']') >> *(omit[Csv::comma_] >> Csv::field_);
        reads_section_heading_.name("read_section_heading_");
        reads_section_ = reads_section_heading_ >> *(section_line_ >> Csv::crlf_);
        reads_section_.name("reads_section_");
        manifests_section_heading_ = char_('[') >> string("Manifests") >> char_(']') >> *(omit[Csv::comma_] >> Csv::field_);
        manifests_section_ = manifests_section_heading_ >> *(section_line_ >> Csv::crlf_);
        manifests_section_.name("manifests_section_");
        settings_section_heading_ = char_('[') >> string("Settings") >> char_(']') >> *(omit[Csv::comma_] >> Csv::field_);
        settings_section_heading_.name("setting_section_header_");
        settings_section_ = settings_section_heading_ >> *(section_line_ >> Csv::crlf_);
        settings_section_.name("settings_section_");

        // unlike the real csv, sample sheets have a special treatment of commented lines
        comment_line_ = char_('#') >> *(char_(0, 0x09) | char_(0x0b, 0x0c) | char_(0x0e, 126)) >> Csv::crlf_;

        name_ = Csv::field_;

        sample_id_  = Csv::field_[_val = bind(&checkIllegalCharacters, "Sample_ID", _1)];
        reference_ = Csv::field_;


        barcode_char_ = char_('A')[_val = 'A'] | char_('a')[_val = 'A'] |
                        char_('C')[_val = 'C'] | char_('c')[_val = 'C'] |
                        char_('G')[_val = 'G'] | char_('g')[_val = 'G'] |
                        char_('T')[_val = 'T'] | char_('t')[_val = 'T'] |
                        char_('N')[_val = 'N'] | char_('n')[_val = 'N'];
        barcode_sequence_ = *barcode_char_;
        description_ = Csv::field_;

        section_line_ = !char_("[") >> Csv::field_ >> *(omit[Csv::comma_] >> Csv::field_);


        data_section_heading_ = char_('[') >> string("Data") >> char_(']') >> *(omit[Csv::comma_] >> Csv::field_);
    }

    // placeholder for operator name extracted at runtime from header section
    std::string operatorName_;
    // placeholder for project name extracted at runtime from header section
    std::string projectName_;

    // sample sheet declarations
    bs::qi::rule<Iterator, bs::unused_type()> comment_line_;
    bs::qi::rule<Iterator, bs::unused_type()> name_;

    bs::qi::rule<Iterator, char()> barcode_char_;
    bs::qi::rule<Iterator, std::string()> barcode_sequence_;

    bs::qi::rule<Iterator, std::string()> project_;
    bs::qi::rule<Iterator, std::string()> operator_;
    bs::qi::rule<Iterator, std::string()> description_;
    bs::qi::rule<Iterator, std::string()> reference_;
    bs::qi::rule<Iterator, std::string()> sample_id_;

    bs::qi::rule<Iterator, bs::unused_type()> section_heading_;
    bs::qi::rule<Iterator, bs::unused_type()> section_line_;

    bs::qi::rule<Iterator, bs::unused_type()> header_section_heading_;
    bs::qi::rule<Iterator, bs::unused_type()> header_section_;

    bs::qi::rule<Iterator, bs::unused_type()> reads_section_heading_;
    bs::qi::rule<Iterator, bs::unused_type()> reads_section_;

    bs::qi::rule<Iterator, bs::unused_type()> manifests_section_heading_;
    bs::qi::rule<Iterator, bs::unused_type()> manifests_section_;

    bs::qi::rule<Iterator, bs::unused_type()> settings_section_heading_;
    bs::qi::rule<Iterator, bs::unused_type()> settings_section_;

    bs::qi::rule<Iterator, bs::unused_type()> data_section_heading_;
};

} // namespace demultiplexing
} // namespace isaac

#endif //iSAAC_DEMULTIPLEXING_MISEQ_SAMPLE_SHEET_CSV_GRAMMAR_HH

