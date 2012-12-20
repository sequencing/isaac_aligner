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
 ** \file ReorderReferenceOptions.cpp
 **
 ** Command line options for 'isaac-reorder-reference'
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_COMMON_REORDER_REFERENCE_OPTIONS_HH
#define iSAAC_COMMON_REORDER_REFERENCE_OPTIONS_HH

#include <string>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace isaac
{
namespace options
{

class ReorderReferenceOptions : public common::Options
{
public:
    boost::filesystem::path sortedReferenceXml_;
    std::string newOrderString_;
    unsigned basesPerLine_;
    std::vector<std::string> newOrder_;
    boost::filesystem::path newXmlPath_;
    boost::filesystem::path newFaPath_;

public:
    ReorderReferenceOptions();

    common::Options::Action parse(int argc, char *argv[]);

private:
    std::string usagePrefix() const {return "isaac-reorder-reference";}
    void postProcess(boost::program_options::variables_map &vm);
    void verifyMandatoryPaths(boost::program_options::variables_map &vm);
};

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_COMMON_REORDER_REFERENCE_OPTIONS_HH
