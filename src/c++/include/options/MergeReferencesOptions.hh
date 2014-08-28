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
 ** \file MergeReferencesOptions.cpp
 **
 ** Command line options for mergeReferences
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_MERGE_REFERENCES_OPTIONS_HH
#define iSAAC_OPTIONS_MERGE_REFERENCES_OPTIONS_HH

#include <string>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace isaac
{
namespace options
{

class MergeReferencesOptions  : public common::Options
{
public:
    std::vector<boost::filesystem::path> filesToMerge_;
    boost::filesystem::path outputFilePath_;

public:
    MergeReferencesOptions();

    common::Options::Action parse(int argc, char *argv[]);

private:
    std::string usagePrefix() const {return "mergeReferences";}
    void postProcess(boost::program_options::variables_map &vm);
    void verifyMandatoryPaths(boost::program_options::variables_map &vm);
};

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_MERGE_REFERENCES_OPTIONS_HH
