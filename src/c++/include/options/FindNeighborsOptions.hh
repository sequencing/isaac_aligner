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
 ** \file FindNeighborsOptions.hh
 **
 ** \brief Command line options for 'findNeighbors'
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_COMMON_FIND_NEIGHBORS_OPTIONS_HH
#define iSAAC_COMMON_FIND_NEIGHBORS_OPTIONS_HH

#include <string>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace isaac
{
namespace options
{

class FindNeighborsOptions : public isaac::common::Options
{
public:
    FindNeighborsOptions();
private:
    std::string usagePrefix() const {return "findNeighbors";}
    void postProcess(boost::program_options::variables_map &vm);
public:
    boost::filesystem::path inputFile;
    boost::filesystem::path outputFile;
    boost::filesystem::path outputDirectory;
    boost::filesystem::path tempFile;
    unsigned jobs;
};

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_COMMON_FIND_NEIGHBORS_OPTIONS_HH
