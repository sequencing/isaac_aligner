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
        allTallies_.back().at(iteration).path_ =
            tempDirectory_ / (boost::format(MATCH_FILE_NAME_TEMPLATE) % flowcellId % lane % tile % iteration).str();
    }
}

size_t MatchTally::getMaxFilePathLength() const
{
    std::size_t ret = 0;
    BOOST_FOREACH(const FileTallyList &list, allTallies_)
    {
        BOOST_FOREACH(const FileTally &tally, list)
        {
            if (ret < tally.path_.string().size())
            {
                ret = tally.path_.string().size();
            }
        }
    }
    return ret;
}

const boost::filesystem::path &MatchTally::getTilePath(const unsigned iteration, const unsigned tileIndex) const
{
    return allTallies_[tileIndex][iteration].path_;
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
    assert(!ft.path_.empty());
    ++ft.matchCount_;
    ++ft.barcodeTally_.at(barcodeIndex);
}

} // namespace alignemnt
} // namespace isaac
