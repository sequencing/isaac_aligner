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
 ** \file FastqFlowcell.hh
 **
 ** Generate flowcell object out of BaseCalls/config.xml
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_ALIGN_OPTIONS_FASTQ_FLOWCELL_HH
#define iSAAC_OPTIONS_ALIGN_OPTIONS_FASTQ_FLOWCELL_HH

#include "flowcell/Layout.hh"
#include "reference/ReferenceMetadata.hh"

namespace isaac
{
namespace options
{
namespace alignOptions
{

struct FastqFlowcellInfo
{
    std::string flowcellId_;
    std::pair<unsigned, unsigned> readLengths_;
    std::vector<unsigned> lanes_;

    const std::vector<unsigned> &getLanes() const
    {
        return lanes_;
    }
};


class FastqFlowcell : boost::noncopyable
{
public:
    static flowcell::Layout createFilteredFlowcell(
        const bool detectSimpleIndels,
        const std::string &tilesFilter,
        const boost::filesystem::path &baseCallsDirectory,
        const bool compressed,
        const unsigned laneNumberMax,
        std::string useBasesMask,
        const bool allowVariableFastqLength,
        const std::string &seedDescriptor,
        const unsigned seedLength,
        const reference::ReferenceMetadataList &referenceMetadataList,
        unsigned &firstPassSeeds);

private:
    struct FastqPathPair
    {
        unsigned lane_;
        boost::filesystem::path r1Path_;
        boost::filesystem::path r2Path_;
    };
    typedef std::vector<FastqPathPair> FastqPathPairList;

    static FastqPathPairList findFastqPathPairs(
        const bool compressed,
        const unsigned laneNumberMax,
        const boost::filesystem::path &baseCallsDirectory);
    static FastqFlowcellInfo parseFastqFlowcellInfo(
        const FastqPathPair &laneFilePaths);
    static FastqFlowcellInfo parseFastqFlowcellInfo(
        const FastqPathPairList &laneFilePaths,
        const bool allowVariableFastqLength);

};

inline std::ostream &operator<< (std::ostream &os, const FastqFlowcellInfo &fcInfo)
{
    os << "FastqFlowcellInfo("<<
        fcInfo.flowcellId_ << "," <<
        fcInfo.readLengths_.first << ":" <<
        fcInfo.readLengths_.second << ",[";

    BOOST_FOREACH(const unsigned lane, fcInfo.getLanes())
    {
        os << lane << " ";
    }
    return os << "])";
}

} // namespace alignOptions
} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_ALIGN_OPTIONS_FASTQ_FLOWCELL_HH
