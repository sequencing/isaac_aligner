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
 ** \file FastqFlowcell.hh
 **
 ** Generate flowcell object out of BaseCalls/laneX.bam files
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_ALIGN_OPTIONS_BAM_FLOWCELL_HH
#define iSAAC_OPTIONS_ALIGN_OPTIONS_BAM_FLOWCELL_HH

#include "flowcell/Layout.hh"
#include "reference/ReferenceMetadata.hh"

namespace isaac
{
namespace options
{
namespace alignOptions
{

struct BamFlowcellInfo
{
    std::string flowcellId_;
    std::pair<unsigned, unsigned> readLengths_;
    std::vector<unsigned> lanes_;

    const std::vector<unsigned> &getLanes() const
    {
        return lanes_;
    }
};


class BamFlowcell : boost::noncopyable
{
    static const unsigned MAX_LANE_NUMBER = 8;

public:
    static flowcell::Layout createFilteredFlowcell(
        const bool detectSimpleIndels,
        const std::string &tilesFilter,
        const boost::filesystem::path &baseCallsDirectory,
        const unsigned laneNumberMax,
        std::string useBasesMask,
        const bool allowVariableReadLength,
        const std::string &seedDescriptor,
        const unsigned seedLength,
        const reference::ReferenceMetadataList &referenceMetadataList,
        unsigned &firstPassSeeds);

private:
    struct BamPath
    {
        unsigned lane_;
        boost::filesystem::path path_;
    };
    static BamPath findBamPath(
        const boost::filesystem::path &baseCallsDirectory);
    static BamFlowcellInfo parseBamFlowcellInfo(
        const BamPath &laneFilePaths,
        const bool allowVariableReadLength,
        const bool allowMixedFlowcells);

};

inline std::ostream &operator<< (std::ostream &os, const BamFlowcellInfo &fcInfo)
{
    os << "BamFlowcellInfo("<<
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

#endif // #ifndef iSAAC_OPTIONS_ALIGN_OPTIONS_BAM_FLOWCELL_HH
