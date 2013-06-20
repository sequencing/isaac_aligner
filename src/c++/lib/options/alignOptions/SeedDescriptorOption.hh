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
 ** \file SeedDescriptorOption.hh
 **
 ** seeds option parsing
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_SEED_DESCRIPTOR_OPTION_HH
#define iSAAC_OPTIONS_SEED_DESCRIPTOR_OPTION_HH

#include "basecalls/ConfigXml.hh"
#include "flowcell/ReadMetadata.hh"

namespace isaac
{
namespace options
{
namespace alignOptions
{

alignment::SeedMetadataList parseSeedDescriptor(
    const bool detectSimpleIndels,
    const std::vector<flowcell::ReadMetadata> &readMetadataList,
    const std::string &seedDescriptor,
    const unsigned seedLength,
    unsigned &firstPassSeeds);

} // namespace alignOptions
} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_SEED_DESCRIPTOR_OPTION_HH
