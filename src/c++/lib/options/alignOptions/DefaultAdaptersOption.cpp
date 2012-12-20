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
 ** \file DefaultAdaptersOption.cpp
 **
 ** default-adapters option parsing
 **
 ** \author Roman Petrovski
 **/
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/mpl/assert.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "flowcell/SequencingAdapterListGrammar.hpp"

#include "DefaultAdaptersOption.hh"


namespace isaac
{
namespace options
{

const flowcell::SequencingAdapterMetadataList parseDefaultAdapters (const std::string &defaultAdapters)
{
    flowcell::SequencingAdapterMetadataList result;
    std::string::const_iterator parseIt(defaultAdapters.begin());
    const std::string::const_iterator parseEnd(defaultAdapters.end());
    flowcell::SequencingAdapterListGrammar<std::string::const_iterator> parser;

    if (!boost::spirit::qi::parse(parseIt, parseEnd, parser, result) || parseEnd != parseIt)
    {
        const boost::format message = boost::format("\n   *** Could not parse the default-adapters '%s' at: %s ***\n") %
            defaultAdapters % defaultAdapters.substr(parseIt - defaultAdapters.begin());
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    std::string adapters;
    BOOST_FOREACH(const flowcell::SequencingAdapterMetadata &adapter, result)
    {
        adapters += (boost::format("%s,") % adapter).str();
    }
    ISAAC_THREAD_CERR << "default-adapters: " << adapters << std::endl;

    return result;
}

} //namespace option
} // namespace isaac
