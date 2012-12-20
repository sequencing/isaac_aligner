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
 ** \file Matchfilter.cpp
 **
 ** Provides a filtering mechanism for the matches, based on the number of
 ** mismatches in each block.
 ** 
 ** \author Come Raczy
 **/

#include <map>
#include <boost/assign.hpp>
#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "alignment/MatchFilter.hh"

namespace isaac
{
namespace alignment
{

std::vector<bool> MatchFilter::getUse(const std::string permutation)
{
    typedef std::map<std::string, std::vector<bool> > UseMap;
    using boost::assign::list_of;
    //                                             (0,0)  (0,1)  (0,2)  (0,3)  (1,0)  (1,1)  (1,2)  (2,0)
    static const std::vector<bool> abcd = list_of( true)( true)( true)(false)( true)( true)(false)( true);
    static const std::vector<bool> bcda = list_of(false)( true)( true)(false)(false)( true)(false)(false);
    static const std::vector<bool> cdab = list_of(false)( true)( true)(false)(false)( true)(false)(false);
    static const std::vector<bool> acbd = list_of(false)(false)(false)(false)(false)( true)(false)(false);
    static const std::vector<bool> bdac = list_of(false)(false)(false)(false)(false)( true)(false)(false);
    static const std::vector<bool> adbc = list_of(false)(false)(false)(false)(false)( true)(false)(false);
    static const UseMap useMap = boost::assign::map_list_of
        ("ABCD", abcd)("BCDA", bcda)("CDAB", cdab)("ACBD", acbd)("BDAC", bdac)("ADBC", adbc);
    const UseMap::const_iterator found = useMap.find(permutation);
    if (useMap.end() == found)
    {
        using boost::format;
        using isaac::common::InvalidParameterException;
        const format message = (format("Unsupported permutation: %s") % permutation);
        BOOST_THROW_EXCEPTION(InvalidParameterException(message.str()));
    }
    return found->second;
}

} // namespace alignment
} // namespace isaac
