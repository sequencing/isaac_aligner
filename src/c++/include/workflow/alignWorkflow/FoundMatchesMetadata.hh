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
 ** \file FoundMatchesMetadata.hh
 **
 ** \brief auxiliary class
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FOUND_MATCHES_METADATA_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_FOUND_MATCHES_METADATA_HH

#include "alignment/MatchTally.hh"
#include "alignment/MatchDistribution.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

struct FoundMatchesMetadata
{
    FoundMatchesMetadata(const boost::filesystem::path &tempDirectory,
                         const flowcell::BarcodeMetadataList &barcodeMetadataList,
                         const unsigned maxIterations,
                         const reference::SortedReferenceXmlList &sortedReferenceXmlList):
                             matchTally_(maxIterations, tempDirectory, barcodeMetadataList),
                             matchDistribution_(sortedReferenceXmlList)
    {

    }

    flowcell::TileMetadataList tileMetadataList_;
    alignment::MatchTally matchTally_;
    alignment::MatchDistribution matchDistribution_;

    void addTile(const flowcell::TileMetadata& tile)
    {
        flowcell::TileMetadata tileWithNewIndex(tile, tileMetadataList_.size());
        tileMetadataList_.push_back(tileWithNewIndex);
        matchTally_.addTile(tileWithNewIndex);
    }
};

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FOUND_MATCHES_METADATA_HH
