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
 ** \file MatchWriter.cpp
 **
 ** Abstract component to write matches.
 **
 ** \author Roman Petrovski
 **/

#include <fstream>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/system/error_code.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "alignment/Match.hh"
#include "io/MatchWriter.hh"

namespace isaac
{
namespace io
{

TileMatchWriter::TileMatchWriter(
    alignment::MatchTally &matchTally,
    const unsigned maxTiles,
    const unsigned maxTileIndex)
    : matchTally_(matchTally),
      tileFileBuffers_(
          maxTiles, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary,
          matchTally_.getMaxFilePathLength()),
      currentIteration_(-1U),
      tileMutexes_(maxTiles)
{
    ISAAC_THREAD_CERR << "Resizing tileStreams to " << maxTileIndex + 1 << std::endl;
    tileStreams_.resize(maxTileIndex + 1);
    ISAAC_THREAD_CERR << "Resized tileStreams to " << tileStreams_.size() << std::endl;

    // create ostream objects for all possible tile indexes
    for (unsigned i = 0; i < tileStreams_.size(); ++i)
    {
        tileStreams_.at(i) = boost::shared_ptr<std::ostream>(new std::ostream(0));
        tileMutexes_.push_back(new boost::mutex);
    }
}

void TileMatchWriter::reopen(const unsigned iteration, const TileMetadataList &tileMetadataList)
{
    ISAAC_ASSERT_MSG(tileFileBuffers_.size() >= tileMetadataList.size(), "Can't be more tiles than initially promised");
    tileFileBuffers_.clear();

    // just in case, disassociate the buffers from ostreams that we will not need
    for (unsigned i = 0; i < tileStreams_.size(); ++i)
    {
        tileStreams_.at(i)->rdbuf(0);
    }

    BOOST_FOREACH(const flowcell::TileMetadata &tile, tileMetadataList)
    {
        const boost::filesystem::path &filePath = matchTally_.getTilePath(iteration, tile.getIndex());

        // this reopens the file handle with new filePath and associates the needed ostream at the tile index position
        tileStreams_.at(tile.getIndex())->rdbuf(tileFileBuffers_.get(filePath, io::FileBufWithReopen::SequentialOnce));
    }
    currentIteration_ = iteration;
}

void TileMatchWriter::write(const SeedId &seedId, const ReferencePosition &referencePosition)
{
    const size_t index = seedId.getTile();

    ISAAC_ASSERT_MSG(0 != tileStreams_.at(index), "Reopen was supposed to create an ostream at this position");
    std::ostream &os = *tileStreams_.at(index);
    const alignment::Match match(seedId, referencePosition);

    boost::lock_guard<boost::mutex> lock(tileMutexes_[index]);
    if (!os.write(reinterpret_cast<const char *>(&match), sizeof(match)))
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to write match into %s") %
            matchTally_.getTilePath(currentIteration_, index)).str()));
    }
    matchTally_(currentIteration_, index, seedId.getBarcode());
}

} //namespace io
} //namespace isaac
