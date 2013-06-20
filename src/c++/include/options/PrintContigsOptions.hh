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
 ** \file SortReferenceOptions.cpp
 **
 ** Command line options for 'sortReference'
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_COMMON_PRINT_CONTIGS_OPTIONS_HH
#define iSAAC_COMMON_PRINT_CONTIGS_OPTIONS_HH

#include <string>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace isaac
{
namespace options
{

class PrintContigsOptions : public isaac::common::Options
{
public:
    PrintContigsOptions();
private:
    std::string usagePrefix() const {return "printContigs";}
    void postProcess(boost::program_options::variables_map &vm);
public:
    boost::filesystem::path originalMetadataPath;
    boost::filesystem::path genomeFile;
};

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_COMMON_PRINT_CONTIGS_OPTIONS_HH
