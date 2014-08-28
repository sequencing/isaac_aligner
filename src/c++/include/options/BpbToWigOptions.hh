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
 ** \file BpbToWigOptions.cpp
 **
 ** Command line options for extractNeighbors
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_BPB_TO_WIG_OPTIONS_HH
#define iSAAC_OPTIONS_BPB_TO_WIG_OPTIONS_HH

#include <string>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace isaac
{
namespace options
{

class BpbToWigOptions  : public common::Options
{
public:
    boost::filesystem::path sortedReferenceMetadata_;
    std::string outputFormatString_;
    boost::filesystem::path inputFilePath_;

public:
    BpbToWigOptions();

    common::Options::Action parse(int argc, char *argv[]);

private:
    std::string usagePrefix() const {return "bpbToWig";}
    void postProcess(boost::program_options::variables_map &vm);
    void verifyMandatoryPaths(boost::program_options::variables_map &vm);
};

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_BPB_TO_WIG_OPTIONS_HH
