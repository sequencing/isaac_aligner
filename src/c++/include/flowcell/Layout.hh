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
 ** \file Layout.hh
 **
 ** Layout of a flowcell.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_FLOWCELL_LAYOUT_HH
#define iSAAC_FLOWCELL_LAYOUT_HH

#include <vector>
#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "alignment/SeedMetadata.hh"

namespace isaac
{
namespace flowcell
{

class Layout
{
public:

    enum Format
    {
        Bam,
        Bcl,
        BclGz,
        Fastq,
        FastqGz
    };

    Layout(const boost::filesystem::path &baseCallsDirectory,
           const Format format,
           const std::vector<unsigned> &barcodeCycles,
           const flowcell::ReadMetadataList &readMetadataList,
           const alignment::SeedMetadataList &seedMetadataList,
           const std::string &flowcellId);

    const boost::filesystem::path &getBaseCallsPath() const {return baseCallsPath_;}
    Format getFormat() const {return format_;}
    const std::string &getFlowcellId() const {return flowcellId_;}
    void setFlowcellId(const std::string &flowcellId) {flowcellId_ = flowcellId;}
    void addTile(const unsigned int lane, const unsigned int tile){
        laneTiles_.resize(std::max<unsigned int>(lane+1, laneTiles_.size()));
        laneTiles_.at(lane).insert(tile);
    }

    bool isEmptyLane(const unsigned int lane) const {
        return laneTiles_.at(lane).empty();
    }

    bool hasLane(const unsigned lane) const
    {
        return !isEmptyLane(lane);
    }

    std::vector<unsigned int> getLaneIds() const
    {
        std::vector<unsigned int> ret;
        std::remove_copy_if(boost::make_counting_iterator(LaneTiles::size_type(0)),
                            boost::make_counting_iterator(laneTiles_.size()),
                            std::back_inserter(ret),
                            boost::bind(&Layout::isEmptyLane, this, _1));
        return ret;
    }

    std::vector<unsigned int> getTileIds(const unsigned int lane) const
    {
        const Tiles &tiles(laneTiles_.at(lane));
        const std::vector<unsigned int> ret(tiles.begin(), tiles.end());
        return ret;
    }

    const std::vector<unsigned> &getBarcodeCycles() const {return barcodeCycles_;}
    unsigned getBarcodeLength() const {return barcodeCycles_.size();}
    const flowcell::ReadMetadataList &getReadMetadataList() const {return readMetadataList_;}
    const alignment::SeedMetadataList &getSeedMetadataList() const {return seedMetadataList_;}
    const std::vector<unsigned> &getAllCycleNumbers() const {return allCycleNumbers_;}

    unsigned getIndex() const {return index_;}
    void setIndex(unsigned index) {index_ = index;}

    const std::pair<std::string, std::string> &getSoftwareVersion() const {return softwareVersion_;}

    void getFiltersFilePath(
        const unsigned tile,
        const unsigned lane,
        boost::filesystem::path &result) const
    {
        ISAAC_ASSERT_MSG(Bcl == format_ || BclGz == format_, "getFiltersFilePath is only allowed for bcl flowcells");
        return getFiltersFilePath(tile, lane, getBaseCallsPath(), result);
    }

    void getFiltersFilePath(
        const unsigned tile,
        const unsigned lane,
        const boost::filesystem::path &baseCallsPath,
        boost::filesystem::path &result) const;

    void getPositionsFilePath(
        const unsigned tile,
        const unsigned lane,
        boost::filesystem::path &result) const
    {
        ISAAC_ASSERT_MSG(Bcl == format_ || BclGz == format_, "getPositionsFilePath is only allowed for bcl flowcells");
        return getPositionsFilePath(tile, lane, getBaseCallsPath(), result);
    }

    void getPositionsFilePath(
        const unsigned tile,
        const unsigned lane,
        const boost::filesystem::path &baseCallsPath,
        boost::filesystem::path &result) const;

    void getBamFilePath(
        boost::filesystem::path &result) const
    {
        ISAAC_ASSERT_MSG(Bam == format_, "getBamFilePath is only allowed for bam flowcells");
        getBamFilePath(getBaseCallsPath(), result);
    }

    static void getBamFilePath(
        const boost::filesystem::path &baseCallsPath,
        boost::filesystem::path &result);

    std::size_t getBamFileSize() const;

    void getFastqFilePath(
        const unsigned read,
        const unsigned lane,
        boost::filesystem::path &result) const
    {
        ISAAC_ASSERT_MSG(Fastq == format_ || FastqGz == format_, "getFastqFilePath is only allowed for fastq flowcells");
        return getFastqFilePath(read, lane, getBaseCallsPath(), getFormat() == flowcell::Layout::FastqGz, result);
    }

    static void getFastqFilePath(
        const unsigned read,
        const unsigned lane,
        const boost::filesystem::path &baseCallsPath,
        const bool compressed,
        boost::filesystem::path &result);

    void getBclFilePath(
        const unsigned tile,
        const unsigned lane,
        const unsigned cycle,
        boost::filesystem::path &result) const
    {
        ISAAC_ASSERT_MSG(Bcl == format_ || BclGz == format_, "getBclFilePath is only allowed for bcl flowcells");
        getBclFilePath(tile, lane, getBaseCallsPath(), cycle, BclGz == getFormat(), result);
    }

    static void getBclFilePath(
        const unsigned tile,
        const unsigned lane,
        const boost::filesystem::path &baseCallsPath,
        const unsigned cycle,
        const bool compressed,
        boost::filesystem::path &result);

    static boost::filesystem::path getLongestFilterFilePath(
        const std::vector<Layout> &flowcellLayoutList)
    {
        boost::filesystem::path ret;
        BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
        {
            boost::filesystem::path filterFile;
            flowcell.getFiltersFilePath(
                            maxTileNumber_, maxLaneNumber_, flowcell.getBaseCallsPath(), filterFile);
            if (ret.string().size() < filterFile.string().size())
            {
                ret = filterFile;
            }
        }
        return ret;
    }

    static boost::filesystem::path getLongestPositionsFilePath(
        const std::vector<Layout> &flowcellLayoutList)
    {
        boost::filesystem::path ret;
        BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
        {
            boost::filesystem::path posFile;
            flowcell.getPositionsFilePath(
                            maxTileNumber_, maxLaneNumber_, flowcell.getBaseCallsPath(), posFile);
            if (ret.string().size() < posFile.string().size())
            {
                ret = posFile;
            }
        }
        return ret;
    }


private:
    static const unsigned maxReadNumber_ = 2;
    static const unsigned maxLaneNumber_ = 8;
    static const unsigned maxTileNumber_ = 9999;
    static const unsigned maxCycleNumber_ = 9999;

    boost::filesystem::path baseCallsPath_;
    Format format_;
    std::vector<unsigned> barcodeCycles_;
    std::string flowcellId_;
    // vector of sets at each lane number position.
    // Lowest lane number to date is 1, so position 0 must have an empty tile set
    typedef std::set<unsigned int> Tiles;
    typedef std::vector<Tiles >LaneTiles;
    LaneTiles laneTiles_;
    flowcell::ReadMetadataList readMetadataList_;
    alignment::SeedMetadataList seedMetadataList_;
    std::vector<unsigned> allCycleNumbers_;
    unsigned index_;
    /// The version of the software used to produce the tile
    std::pair<std::string, std::string> softwareVersion_;
    std::pair<unsigned, unsigned> softwareMajorMinor_;

};

typedef std::vector<flowcell::Layout> FlowcellLayoutList;

inline unsigned getMaxReadLength(const flowcell::FlowcellLayoutList &flowcellLayoutList)
{
    unsigned ret = 0;
    BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
    {
        ret = std::max(ret, getMaxReadLength(flowcell.getReadMetadataList()));
    }
    return ret;
}

inline unsigned getMaxReadCount(const FlowcellLayoutList &flowcellLayoutList)
{
    return std::max_element(
        flowcellLayoutList.begin(), flowcellLayoutList.end(),
        boost::bind(&ReadMetadataList::size, boost::bind(&Layout::getReadMetadataList, _1)))->getReadMetadataList().size();
}

inline unsigned getMaxReadLength(const flowcell::FlowcellLayoutList &flowcellLayoutList, const unsigned readIndex)
{
    unsigned ret = 0;
    BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
    {
        ret = std::max(ret,
            flowcell.getReadMetadataList().size() > readIndex ? flowcell.getReadMetadataList()[readIndex].getLength() : 0);
    }
    return ret;
}

inline size_t getMaxSeedsPerRead(const flowcell::FlowcellLayoutList &flowcellLayoutList)
{
    size_t ret = 0;
    BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
    {
        BOOST_FOREACH(const flowcell::ReadMetadata &read, flowcell.getReadMetadataList())
        {
            ret = std::max<size_t>(ret, std::count_if(flowcell.getSeedMetadataList().begin(),
                                                      flowcell.getSeedMetadataList().end(),
                                                      boost::bind(&alignment::SeedMetadata::getReadIndex, _1) ==
                                                          read.getIndex()));
        }
    }
    return ret;
}

inline unsigned getMaxTotalReadLength(const FlowcellLayoutList &flowcellLayoutList)
{
    unsigned ret = 0;
    BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
    {
        ret = std::max(ret, getTotalReadLength(flowcell.getReadMetadataList()));
    }
    return ret;
}

inline unsigned getMinTotalReadLength(const FlowcellLayoutList &flowcellLayoutList)
{
    unsigned ret = std::numeric_limits<unsigned>::max();
    BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
    {
        ret = std::min(ret, getTotalReadLength(flowcell.getReadMetadataList()));
    }
    return ret;
}

inline unsigned getMaxBarcodeLength(const FlowcellLayoutList &flowcellLayoutList)
{
    unsigned ret = 0;
    BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
    {
        ret = std::max(ret, flowcell.getBarcodeLength());
    }
    return ret;
}

inline boost::filesystem::path getLongestBaseCallsPath(const flowcell::FlowcellLayoutList &flowcellLayoutList)
{
    boost::filesystem::path ret;
    BOOST_FOREACH(const flowcell::Layout &flowcellLayout, flowcellLayoutList)
    {
        if (ret.string().size() < flowcellLayout.getBaseCallsPath().string().size())
        {
            ret = flowcellLayout.getBaseCallsPath();
        }
    }
    return ret;
}

inline std::ostream & operator << (std::ostream &os, const Layout &layout)
{
    os << "Layout("
              << layout.getFlowcellId() << ", "
              << layout.getSoftwareVersion().first << "-"
              << layout.getSoftwareVersion().second << "[";
    BOOST_FOREACH(const unsigned cycle, layout.getAllCycleNumbers())
    {
        os << "," << cycle;
    }
    return os
              << "])";
}

} // namespace flowcell
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_LAYOUT_HH
