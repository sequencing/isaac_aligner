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

#ifndef iSAAC_IO_MATCH_WRITER_HH
#define iSAAC_IO_MATCH_WRITER_HH

#include <iostream>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "alignment/SeedId.hh"
#include "alignment/MatchTally.hh"
#include "flowcell/TileMetadata.hh"
#include "io/FileBufCache.hh"
#include "reference/ReferenceKmer.hh"

namespace isaac
{
namespace io
{

/**
 ** \brief a component that encapsulates the binning of matches into several
 ** files.
 **
 ** In this implementation, the binning is done per tile and
 ** per mask and per iteration. This is an implicit coupling to the structure of the MatchFinder
 ** workflow. MatchWriter holds a separate stream for each tile. On each write the referenced
 ** matchTally is updated.
 **/
class TileMatchWriter: boost::noncopyable, boost::ptr_vector<io::FileBufWithReopen>
{
    typedef alignment::SeedId SeedId;
    typedef isaac::reference::ReferencePosition ReferencePosition;
public:
    typedef flowcell::TileMetadata TileMetadata;
    typedef std::vector<TileMetadata> TileMetadataList;
    TileMatchWriter(
        alignment::MatchTally &matchTally,
        const unsigned maxTiles,
        const unsigned maxTileIndex);

    /**
     * \brief Switches to a new set of tile files based on the iteration supplied
     */
    void reopen(const unsigned iteration,
                const TileMetadataList &tileMetadataList);

    /**
     * \brief Extracts the tile id from the SeedId and write the match to the
     *  appropriate stream
     */
    void write(const SeedId &seedId, const ReferencePosition &referencePosition);

private:
    alignment::MatchTally &matchTally_;
    io::FileBufCache<io::FileBufWithReopen> tileFileBuffers_;
    std::vector<boost::shared_ptr<std::ostream> > tileStreams_;
    unsigned currentIteration_;
    boost::ptr_vector<boost::mutex> tileMutexes_;
};
} //namespace io
} //namespace isaac

#endif // #ifndef iSAAC_IO_MATCH_WRITER_HH
