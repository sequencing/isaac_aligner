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
 ** \file ExtractNeighborsOptions.cpp
 **
 ** Command line options for 'isaac-reorder-reference'
 **
 ** \author Roman Petrovski
 **/

#include <algorithm>
#include <string>
#include <vector>
#include <ostream>
#include <fstream>

#include <boost/algorithm/string/regex.hpp>
#include <boost/assign.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lexical_cast.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "options/ExtractNeighborsOptions.hh"



namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;
using common::InvalidOptionException;

ExtractNeighborsOptions::ExtractNeighborsOptions()
{
    namedOptions_.add_options()
        ("reference-genome,r"       , bpo::value<bfs::path>(&sortedReferenceXml_),
                "Full path to the reference genome XML descriptor."
            )
        ("output-file,o"       , bpo::value<bfs::path>(&outputFilePath_),
                "Path for the output file."
            );
}

common::Options::Action ExtractNeighborsOptions::parse(int argc, char *argv[])
{
    const std::vector<std::string> allOptions(argv, argv + argc);
    common::Options::Action ret = common::Options::parse(argc, argv);
    if (RUN == ret)
    {
        ISAAC_THREAD_CERR << "argc: " << argc << " argv: " << boost::join(allOptions, " ") << std::endl;
    }
    return ret;
}

void ExtractNeighborsOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help") ||  vm.count("version"))
    {
        return;
    }

    const std::vector<std::string> requiredOptions = boost::assign::list_of("output-file")("reference-genome");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const boost::format message = boost::format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }

    outputFilePath_ = boost::filesystem::absolute(outputFilePath_);// / "GenomeNeighbors.dat";

    if (boost::filesystem::exists(outputFilePath_))
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException((boost::format("\n   *** Output file already exists. : %s ***\n") % outputFilePath_.string()).str()));
    }
}




} //namespace options
} // namespace isaac
