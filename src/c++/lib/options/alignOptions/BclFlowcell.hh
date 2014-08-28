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
 ** \file BclFlowcell.hh
 **
 ** Generate flowcell object out of BaseCalls/config.xml
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_ALIGN_OPTIONS_BCL_FLOWCELL_HH
#define iSAAC_OPTIONS_ALIGN_OPTIONS_BCL_FLOWCELL_HH

#include "flowcell/Layout.hh"
#include "reference/ReferenceMetadata.hh"

namespace isaac
{
namespace options
{
namespace alignOptions
{

class BclFlowcell : boost::noncopyable
{
public:
    static flowcell::Layout createFilteredFlowcell(
        const bool detectSimpleIndels,
        const std::string &tilesFilter,
        const boost::filesystem::path &baseCallsPath,
        const flowcell::Layout::Format format,
        const bool compressed,
        const unsigned laneNumberMax,
        std::string useBasesMask,
        const std::string &seedDescriptor,
        const unsigned seedLength,
        const reference::ReferenceMetadataList &referenceMetadataList,
        unsigned &firstPassSeeds);

private:

};

} // namespace alignOptions
} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_ALIGN_OPTIONS_BCL_FLOWCELL_HH
