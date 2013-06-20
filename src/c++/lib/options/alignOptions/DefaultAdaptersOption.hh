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
 ** \file DefaultAdaptersOption.hh
 **
 ** default-adapters option parsing
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_DEFAULT_ADAPTERS_OPTION_HH
#define iSAAC_OPTIONS_DEFAULT_ADAPTERS_OPTION_HH

#include "basecalls/ConfigXml.hh"
#include "flowcell/ReadMetadata.hh"

namespace isaac
{
namespace options
{

const flowcell::SequencingAdapterMetadataList parseDefaultAdapters (const std::string &defaultAdapters);

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_DEFAULT_ADAPTERS_OPTION_HH
