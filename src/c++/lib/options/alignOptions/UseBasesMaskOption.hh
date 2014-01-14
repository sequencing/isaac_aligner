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
