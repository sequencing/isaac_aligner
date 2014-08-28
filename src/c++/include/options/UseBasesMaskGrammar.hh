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
 ** \file UseBasesMaskGrammar.hh
 **
 ** boost::spirit grammar for use bases mask parsing
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_USE_BASES_MASK_GRAMMAR_HH
#define iSAAC_OPTIONS_USE_BASES_MASK_GRAMMAR_HH


#include <boost/foreach.hpp>
#include <boost/lambda/if.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_attr_cast.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix1_statements.hpp>
#include <boost/spirit/home/phoenix/statement/if.hpp>
#include <boost/spirit/home/phoenix/statement/throw.hpp>
#include <boost/spirit/home/phoenix/object/static_cast.hpp>
#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/spirit/home/phoenix/bind.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/support/unused.hpp>

namespace isaac
{
namespace options
{

namespace bs=boost::spirit;



template <typename Iterator>
struct UseBasesMaskGrammar : bs::qi::grammar<Iterator, std::vector<std::string >()>
{
    UseBasesMaskGrammar(const std::vector<unsigned int> &readLengths) :
        UseBasesMaskGrammar::base_type(start_), readLengths_(readLengths), currentRead_(0)
    {
        // if there is insufficient number of readLenghts, use the last element to assume 0 for the missing ones
        readLengths_.push_back(0);

        using bs::ascii::char_;
        using bs::qi::_1;
        using bs::qi::_val;
        using bs::uint_;
        using boost::phoenix::at;
        using boost::phoenix::back;
        using boost::phoenix::begin;
        using boost::phoenix::bind;
        using boost::phoenix::construct;
        using boost::phoenix::end;
        using boost::phoenix::if_;
        using boost::phoenix::insert;
        using boost::phoenix::push_back;
        using boost::phoenix::ref;
        using boost::phoenix::size;


        // valid mask characters
        valid_chars_ = char_('I')[_val = 'i'] | char_('i')[_val = 'i']|char_('Y')[_val = 'y'] | char_('y')[_val = 'y'] | char_('N')[_val = 'n'] | char_('n')[_val = 'n'];

        // for some reason (uint_ - uint(0)) results in the attached semantic action not being executed...
        repeat_count_ = uint_[&failIf0][_val = _1];

        // expansion for simple things like y or y10
        read_mask_no_wc_ =
            *(
                valid_chars_[push_back(_val, _1)] >> -repeat_count_[insert(_val, end(_val), _1-1, back(_val))]
            );

        // allow only one star per read mask. Then consume the rest of the read by the read_mask_no_wc
        read_mask_wc_ =
            *(
                valid_chars_[push_back(_val, _1)] >>
                -(
                    repeat_count_[insert(_val, end(_val), _1-1, back(_val))] |
                    (char_('*') >>
                        read_mask_no_wc_
                        [
                             if_(at(ref(readLengths_), ref(currentRead_)) > size(_val) + size(_1))
                                 [insert(_val, end(_val),
                                         at(ref(readLengths_), ref(currentRead_)) - size(_val) - size(_1),
                                         back(_val))]
                        ]
                        [insert(_val, end(_val), begin(_1), end(_1))]
//                                [bind(&appendBoth, _val
//                                            ,
//                                            at(readLengths_, val(currentRead_))-static_cast_<int>(size(_val)+size(_1)),
//                                            back(_val),
//                                            _1
//                                            )]
                    )
                )
            );
        use_bases_mask_ = read_mask_wc_
            [
             //stay on last read length (of 0), dont't go over the vector
             if_((readLengths_.size() - 1) != ref(currentRead_))
                 [ref(currentRead_)++]
            ]
            [push_back(_val, construct<std::string>(begin(_1),end(_1)))] % ',';

        start_ = use_bases_mask_;
    }

    //    static void appendBoth(std::vector<char> &where, int howMany, char what, const std::vector<char>& rest)
    //    {
    //        std::cerr << "where:" << std::string(where.begin(), where.end()) <<
    //                " howMany:" << howMany << " what:" << what << " rest:" << std::string(rest.begin(), rest.end());
    //        where.insert(where.end(), howMany, what);
    //        where.insert(where.end(), rest.begin(), rest.end());
    //    }

    static void failIf0(unsigned int attribute, const bs::unused_type& it, bool &mFlag)
    {
        if (!attribute)
        {
            mFlag = false;
        }
    }
    std::vector<unsigned int> readLengths_;
    unsigned currentRead_;

    bs::qi::rule<Iterator, char()> valid_chars_;
    bs::qi::rule<Iterator, unsigned()> repeat_count_;
    bs::qi::rule<Iterator, std::vector<char>() > read_mask_wc_;
    bs::qi::rule<Iterator, std::vector<char>()> read_mask_no_wc_;
    bs::qi::rule<Iterator, std::vector<std::string >()> use_bases_mask_;
    bs::qi::rule<Iterator, std::vector<std::string >()> start_;
};

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_USE_BASES_MASK_GRAMMAR_HH
