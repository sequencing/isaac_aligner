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
    static const unsigned MAX_LANE_NUMBER = 8;

public:
    static flowcell::Layout createFilteredFlowcell(
        const std::string &tilesFilter,
        const boost::filesystem::path &baseCallsDirectory,
        const flowcell::Layout::Format format,
        std::string useBasesMask,
        const std::string &seedDescriptor,
        const reference::ReferenceMetadataList &referenceMetadataList);

private:
    struct FastqPathPair
    {
        unsigned lane_;
        boost::filesystem::path r1Path_;
        boost::filesystem::path r2Path_;
    };
    typedef std::vector<FastqPathPair> FastqPathPairList;

    static FastqPathPairList findFastqPathPairs(
        const flowcell::Layout::Format format,
        const boost::filesystem::path &baseCallsDirectory);
    static FastqFlowcellInfo parseFastqFlowcellInfo(
        const FastqPathPair &laneFilePaths);
    static FastqFlowcellInfo parseFastqFlowcellInfo(
        const FastqPathPairList &laneFilePaths);

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
