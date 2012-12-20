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
 ** \file MatchTally.hh
 **
 ** Keeps track of the match count in each of the match files produced.
 **
 ** \author Come Raczy
 **/

#include <boost/foreach.hpp>

#include "alignment/MatchTally.hh"

namespace isaac
{
namespace alignment
{

static const char * MATCH_FILE_NAME_TEMPLATE = "%s_s_%d_%04d_%d_matches.dat";

MatchTally::MatchTally(const unsigned maxIterations,
                       const boost::filesystem::path &tempDirectory,
                       const flowcell::BarcodeMetadataList &barcodeMetadataList)
    : maxIterations_(maxIterations),
      barcodeMetadataList_(barcodeMetadataList),
      tempDirectory_(tempDirectory)
{
}

void MatchTally::addTile(const flowcell::TileMetadata& tileMetadata)
{
    allTallies_.push_back(FileTallyList(maxIterations_, FileTally(barcodeMetadataList_.size())));
    for (unsigned iteration = 0; iteration < maxIterations_; ++iteration)
    {
        const std::string &flowcellId = tileMetadata.getFlowcellId();
        const unsigned int lane = tileMetadata.getLane();
        const unsigned int tile = tileMetadata.getTile();
        allTallies_.back().at(iteration).first =
            tempDirectory_ / (boost::format(MATCH_FILE_NAME_TEMPLATE) % flowcellId % lane % tile % iteration).str();
    }
}

size_t MatchTally::getMaxFilePathLength(const boost::filesystem::path &tempDirectory)
{
    // use some ridiculously long flowcell id and other components
    return (tempDirectory / (boost::format(MATCH_FILE_NAME_TEMPLATE) % std::string(256, 'x') % 123 % 1234 % 100).str()).string().size();
}

const boost::filesystem::path &MatchTally::getTilePath(const unsigned iteration, const unsigned tileIndex) const
{
    return allTallies_[tileIndex][iteration].first;
}

const MatchTally::FileTallyList &MatchTally::getFileTallyList(const flowcell::TileMetadata &tileMetadata) const
{
    assert(tileMetadata.getIndex() < allTallies_.size());
    return allTallies_[tileMetadata.getIndex()];
}

/**
 * \brief Updates tile file mapping and match count statistics. Threads are expected to
 *        not update the same tile simultaneously
 */
void MatchTally::operator()(const unsigned iteration, const unsigned tileIndex, const unsigned barcodeIndex)
{
    FileTally &ft = allTallies_.at(tileIndex).at(iteration);
    assert(!ft.first.empty());
    ++ft.second;
    ++ft.barcodeTally_.at(barcodeIndex);
}

} // namespace alignemnt
} // namespace isaac
