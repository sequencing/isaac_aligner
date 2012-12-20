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
 ** \file .cpp
 **
 ** Command line options for extractNeighbors
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_EXTRACT_NEIGHBORS_OPTIONS_HH
#define iSAAC_OPTIONS_EXTRACT_NEIGHBORS_OPTIONS_HH

#include <string>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace isaac
{
namespace options
{

class ExtractNeighborsOptions  : public common::Options
{
public:
    boost::filesystem::path sortedReferenceXml_;
    boost::filesystem::path outputFilePath_;

public:
    ExtractNeighborsOptions();

    common::Options::Action parse(int argc, char *argv[]);

private:
    std::string usagePrefix() const {return "extractNeighbors";}
    void postProcess(boost::program_options::variables_map &vm);
    void verifyMandatoryPaths(boost::program_options::variables_map &vm);
};

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_EXTRACT_NEIGHBORS_OPTIONS_HH
