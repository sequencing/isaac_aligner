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
 ** \file UseBasesMaskOption.hh
 **
 ** use-bases-mask option parsing
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_USE_BASES_MASK_OPTION_HH
#define iSAAC_OPTIONS_USE_BASES_MASK_OPTION_HH

#include "rta/ConfigXml.hh"
#include "flowcell/ReadMetadata.hh"

namespace isaac
{
namespace options
{

struct ParsedUseBasesMask
{
    flowcell::ReadMetadataList dataReads_;
    flowcell::ReadMetadataList indexReads_;
};

ParsedUseBasesMask parseUseBasesMask (
    const std::vector<unsigned int> &readFirstCycles,
    const std::vector<unsigned int> &readLengths,
    const unsigned seedLength,
    const std::string &useBasesMask,
    const boost::filesystem::path &baseCallsDirectory);

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_USE_BASES_MASK_GRAMMAR_HH
