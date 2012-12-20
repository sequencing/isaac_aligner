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
 ** \brief Tracking of the match count in each of the match files produced.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_TALLY_HH
#define iSAAC_ALIGNMENT_MATCH_TALLY_HH

#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>

#include "common/Debug.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "reference/SortedReferenceXml.hh"
#include "oligo/Mask.hh"

namespace isaac
{
namespace alignment
{

namespace bfs = boost::filesystem;


struct MatchTally
{
    struct FileTally
    {
        FileTally(const size_t barcodes = 0) : second(0), barcodeTally_(barcodes){
//            ISAAC_THREAD_CERR << "constructed FileTally for " << barcodeTally_.size() << "barcodes\n";
        }
        unsigned long getBarcodeMatchCount(const unsigned barcode) const
        {
            return barcodeTally_.at(barcode);
        }
        bfs::path first;
        unsigned long second;
        std::vector<unsigned long> barcodeTally_;
    };
    typedef std::vector<FileTally> FileTallyList;

public:
    MatchTally(const unsigned maxIterations,
               const boost::filesystem::path &tempDirectory,
//               const std::vector<flowcell::TileMetadata> &tileMetadataList,
               const flowcell::BarcodeMetadataList &barcodeMetadataList);

    static size_t getMaxFilePathLength(const boost::filesystem::path &tempDirectory);

    const boost::filesystem::path &getTilePath(const unsigned iteration, const unsigned tileIndex) const;

    /// record the match count for each file produced by the matchWriter
    void operator()(const unsigned iteration, const unsigned tileIndex, const unsigned barcodeIndex);

    /// return all the tally for all match files for the given tile
    const FileTallyList &getFileTallyList(const flowcell::TileMetadata &tileMetadata) const;

    MatchTally &operator = (MatchTally that)
    {
        swap(*this, that);
        return *this;
    }

    void addTile(const flowcell::TileMetadata& tile);

private:
    const unsigned maxIterations_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const boost::filesystem::path &tempDirectory_;
    std::vector<FileTallyList> allTallies_;

    friend void swap(MatchTally &one, MatchTally &another);

    template<class Archive> friend void serialize(Archive & ar, MatchTally &, const unsigned int file_version);
};

inline void swap(MatchTally &one, MatchTally &another)
{
    one.allTallies_.swap(another.allTallies_);
}

} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_TALLY_HH
