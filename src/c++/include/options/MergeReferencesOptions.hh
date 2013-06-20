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
