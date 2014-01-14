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
