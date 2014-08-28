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
 ** \file SequencingAdapterListGrammar.hpp
 **
 ** grammar for parsing list of sequencing adapter definitions
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_FLOWCELL_SEQUENCING_ADAPTER_LIST_GRAMMAR_HH
#define iSAAC_FLOWCELL_SEQUENCING_ADAPTER_LIST_GRAMMAR_HH

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_attr_cast.hpp>
//#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
//#include <boost/spirit/include/phoenix1_statements.hpp>
//#include <boost/spirit/include/support_istream_iterator.hpp>
//#include <boost/spirit/home/phoenix/statement/if.hpp>
//#include <boost/spirit/home/phoenix/statement/throw.hpp>
//#include <boost/spirit/home/phoenix/object/static_cast.hpp>
//#include <boost/spirit/home/phoenix/function/function.hpp>
//#include <boost/spirit/home/phoenix/bind.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/support/unused.hpp>

#include "flowcell/SequencingAdapterMetadata.hh"

namespace isaac
{
namespace flowcell
{


namespace bs=boost::spirit;


template <typename Iterator>
struct SequencingAdapterListGrammar : bs::qi::grammar<Iterator, flowcell::SequencingAdapterMetadataList() >
{
    SequencingAdapterListGrammar() :
        SequencingAdapterListGrammar::base_type(start_)
    {

        using bs::ascii::char_;
        using bs::ascii::string;
        //using bs::uint_;
        //using bs::double_;
        using bs::qi::_1;
//        using bs::qi::_2;
//        using bs::qi::_3;
//        using bs::qi::_4;
        using bs::qi::omit;
        using bs::qi::_val;

//        using boost::phoenix::at;
//        using boost::phoenix::back;
//        using boost::phoenix::begin;
//        using boost::phoenix::bind;
        using boost::phoenix::construct;
//        using boost::phoenix::end;
//        using boost::phoenix::if_;
//        using boost::phoenix::insert;
        using boost::phoenix::push_back;
//        using boost::phoenix::ref;
        using boost::phoenix::size;


        star_ = char_('*');
        comma_ = char_(',');

        adapter_char_ = char_('A')[_val = 'A'] | char_('a')[_val = 'A'] |
                        char_('C')[_val = 'C'] | char_('c')[_val = 'C'] |
                        char_('G')[_val = 'G'] | char_('g')[_val = 'G'] |
                        char_('T')[_val = 'T'] | char_('t')[_val = 'T'];
        adapter_sequence_ = adapter_char_ >> adapter_char_ >> adapter_char_ >> adapter_char_ >> adapter_char_ >> *adapter_char_;

        simple_adapter_ =
            adapter_sequence_[_val = construct<flowcell::SequencingAdapterMetadata>(_1, false, size(_1))];
        forward_unbounded_adapter_ =
            adapter_sequence_[_val = construct<flowcell::SequencingAdapterMetadata>(_1, false, 0)] >> star_;
        reverse_unbounded_adapter_ =
            star_ >> adapter_sequence_[_val = construct<flowcell::SequencingAdapterMetadata>(_1, true, 0)];
        adapter_ = forward_unbounded_adapter_ | simple_adapter_ | reverse_unbounded_adapter_ ;


        adapter_list_ = *(adapter_ >> -omit[comma_]);

        standard_macro_ = string("Standard")[_val = STANDARD_ADAPTERS];
        nextera_mp_macro_ = string("NexteraMp")[_val = NEXTERA_MATEPAIR_ADAPTERS];
        nextera_macro_ = string("Nextera")[_val = NEXTERA_STANDARD_ADAPTERS];
        macro_ = standard_macro_ | nextera_mp_macro_ | nextera_macro_;
        start_ = macro_ | adapter_list_;
    }

    bs::qi::rule<Iterator, char()> star_;
    bs::qi::rule<Iterator, char()> comma_;

    bs::qi::rule<Iterator, char()> adapter_char_;
    bs::qi::rule<Iterator, std::string()> adapter_sequence_;

    bs::qi::rule<Iterator, flowcell::SequencingAdapterMetadataList()> standard_macro_;
    bs::qi::rule<Iterator, flowcell::SequencingAdapterMetadataList()> nextera_macro_;
    bs::qi::rule<Iterator, flowcell::SequencingAdapterMetadataList()> nextera_mp_macro_;
    bs::qi::rule<Iterator, flowcell::SequencingAdapterMetadataList()> macro_;
    bs::qi::rule<Iterator, flowcell::SequencingAdapterMetadata()> simple_adapter_;
    bs::qi::rule<Iterator, flowcell::SequencingAdapterMetadata()> forward_unbounded_adapter_;
    bs::qi::rule<Iterator, flowcell::SequencingAdapterMetadata()> reverse_unbounded_adapter_;
    bs::qi::rule<Iterator, flowcell::SequencingAdapterMetadata()> adapter_;
    bs::qi::rule<Iterator, flowcell::SequencingAdapterMetadataList() > adapter_list_;
    bs::qi::rule<Iterator, flowcell::SequencingAdapterMetadataList() > start_;
};

} // namespace flowcell
} // namespace isaac

#endif //iSAAC_FLOWCELL_SEQUENCING_ADAPTER_LIST_GRAMMAR_HH
