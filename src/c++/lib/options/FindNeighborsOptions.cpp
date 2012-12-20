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
 ** \file FindNeighborsOptions.cpp
 **
 ** \brief See FindNeighborsOptions.hh.
 **
 ** \author Come Raczy
 **/

#include <string>
#include <vector>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>

#include "options/FindNeighborsOptions.hh"
#include "common/Exceptions.hh"

namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

FindNeighborsOptions::FindNeighborsOptions()
    : inputFile("")
    , outputDirectory("./")
    , tempFile("Temp/neighbors.dat")
    , jobs(boost::thread::hardware_concurrency())
{
    namedOptions_.add_options()
        ("input-file,i",  bpo::value<bfs::path>(&inputFile),
                          "The input 'SortedReference.xml' file")
        ("output-file,o",  bpo::value<bfs::path>(&outputFile),
                          "The output 'SortedReference.xml' file")
        ("output-directory",  bpo::value<bfs::path>(&outputDirectory),
                          "The location for annotated data files")
        ("temp-file,t", bpo::value<bfs::path>(&tempFile)->default_value(tempFile),
                          "The file where all the kmers with a neighborhood will be written")
        ("jobs,j", bpo::value<unsigned>(&jobs)->default_value(jobs),
                          "Maximum number of compute threads to run in parallel. Parallel sorting will use all cores regardless.")
        ;
}

void FindNeighborsOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help"))
    {
        return;
    }
    using isaac::common::InvalidOptionException;
    using boost::format;
    const std::vector<std::string> requiredOptions = boost::assign::list_of("input-file");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const format message = format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
    typedef std::pair<bfs::path *, std::string> PathOption;
    const std::vector<PathOption> pathOptions = boost::assign::list_of
        (PathOption(&inputFile, "input-file"))
        (PathOption(&tempFile,"output-file"))
        ;
    BOOST_FOREACH(const PathOption &pathOption, pathOptions)
    {
        if(pathOption.first->empty())
        {
            const format message = format("\n   *** The '%s' can't be empty (use '.' for current directory) ***\n") % pathOption.second;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
    const std::vector<PathOption> existingPaths = boost::assign::list_of
        (PathOption(&inputFile, "input-file"))
        ;
    BOOST_FOREACH(const PathOption &pathOption, existingPaths)
    {
        if(!exists(*pathOption.first))
        {
            const format message = format("\n   *** The '%s' does not exist: %s ***\n") % pathOption.second % *pathOption.first;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
}

} // namespace options
} // namespace isaac
