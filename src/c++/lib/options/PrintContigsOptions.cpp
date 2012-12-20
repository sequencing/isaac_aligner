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
 ** \file PrintContigsOptions.cpp
 **
 ** Command line options for 'printContigs'
 **
 ** \author Roman Petrovski
 **/

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include "options/PrintContigsOptions.hh"

namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;

PrintContigsOptions::PrintContigsOptions()
{
    namedOptions_.add_options()
        ("original-metadata"       , bpo::value<boost::filesystem::path>(&originalMetadataPath),
                "If supplied, the additional contig metadata fields will be copied over from it."
            )
        ("genome-file,g",       bpo::value<boost::filesystem::path>(&genomeFile),
                                "Name of the reference genome")
        ;
}

void PrintContigsOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help"))
    {
        return;
    }
    using isaac::common::InvalidOptionException;
    using boost::format;
    const std::vector<std::string> requiredOptions = boost::assign::list_of("genome-file");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const format message = format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
}

} //namespace option
} // namespace isaac
