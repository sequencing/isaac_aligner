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
 ** \file CsvGrammar.hh
 **
 ** SampleSheet.csv grammar definition
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_COMMON_CSV_GRAMMAR_HH
#define iSAAC_COMMON_CSV_GRAMMAR_HH

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

namespace isaac
{
namespace common
{


namespace bs=boost::spirit;


template <typename Iterator>
struct CsvGrammar
{
    CsvGrammar()
    {
        using bs::ascii::char_;
        using bs::qi::_1;
        using boost::phoenix::push_back;
        using bs::qi::_val;

        /**
         * CSV ABNF mainly taken from http://www.ietf.org/rfc/rfc4180.txt with some modification to make sample
         * sheet users life easier
         */
        cr_ = char_('\xd');
        comma_ = char_(',');
        comma_.name(",");
        dquote_ = char_('"');
        lf_ = char_('\xa');

        // unlike the real csv, sample sheets often get mixed line endings
        crlf_ = +(cr_ | lf_);
        textdata_ = char_('\x20', '\x21') | char_('\x23', '\x2b') | char_('\x2d', '\x7e');
        escaped_ = dquote_ >>
            *(textdata_[push_back(_val, _1)]|
              comma_[push_back(_val, _1)]|
              cr_[push_back(_val, _1)]|
              lf_[push_back(_val, _1)]|
              (dquote_>>dquote_)[push_back(_val, '"')]) >>
            dquote_;
        non_escaped_ = *textdata_;
        field_ = escaped_ | non_escaped_;
    }

    // standard CSV declarations
    bs::qi::rule<Iterator, std::string()> field_;
    bs::qi::rule<Iterator, std::string()> escaped_;
    bs::qi::rule<Iterator, std::string()> non_escaped_;
    bs::qi::rule<Iterator, char()> comma_;
    bs::qi::rule<Iterator, char()> cr_;
    bs::qi::rule<Iterator, char()> dquote_;
    bs::qi::rule<Iterator, char()> lf_;
    bs::qi::rule<Iterator, bs::unused_type()> crlf_;
    bs::qi::rule<Iterator, char()> textdata_;
};

} // namespace demultiplexing
} // namespace isaac

#endif //iSAAC_COMMON_CSV_GRAMMAR_HH
