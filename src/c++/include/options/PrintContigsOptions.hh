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
