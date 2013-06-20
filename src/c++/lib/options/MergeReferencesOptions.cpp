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
 ** \file MergeReferencesOptions.cpp
 **
 ** Command line options for 'mergeReferences'
 **
 ** \author Roman Petrovski
 **/

#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "options/MergeReferencesOptions.hh"



namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;
using common::InvalidOptionException;

MergeReferencesOptions::MergeReferencesOptions()
{
    namedOptions_.add_options()
        ("input-file,i"    , bpo::value<std::vector<bfs::path> >(&filesToMerge_),
                "Paths of the files to be merged."
            )
        ("output-file,o"       , bpo::value<bfs::path>(&outputFilePath_), "Path for the output file."
            );
}

common::Options::Action MergeReferencesOptions::parse(int argc, char *argv[])
{
    const std::vector<std::string> allOptions(argv, argv + argc);
    common::Options::Action ret = common::Options::parse(argc, argv);
    if (RUN == ret)
    {
        ISAAC_THREAD_CERR << "argc: " << argc << " argv: " << boost::join(allOptions, " ") << std::endl;
    }
    return ret;
}

void MergeReferencesOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help") ||  vm.count("version"))
    {
        return;
    }

    const std::vector<std::string> requiredOptions = boost::assign::list_of("input-file")("output-file");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const boost::format message = boost::format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
}




} //namespace options
} // namespace isaac
