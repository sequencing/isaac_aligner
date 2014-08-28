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
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/variant.hpp>

#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "alignment/SeedMetadata.hh"

namespace isaac
{
namespace flowcell
{

struct BclFlowcellData
{
    BclFlowcellData() : compressed_(false), patternedFlowcell_(false), tilesPerLaneMax_(0){}
    std::pair<std::string, std::string> softwareVersion_;
    std::pair<unsigned, unsigned> softwareMajorMinor_;
    bool compressed_;
    bool patternedFlowcell_;
    // number of tiles a lane can possibly have regardless of tile filtering
    unsigned tilesPerLaneMax_;
};

struct FastqFlowcellData
{
    FastqFlowcellData(bool compressed) : compressed_(compressed){}
    bool compressed_;
};

struct BamFlowcellData
{

};


class Layout
{
public:
    enum Format
    {
        Bam,
        Bcl,
        BclBgzf,
        Fastq,
    };
    typedef boost::variant<BclFlowcellData, FastqFlowcellData, BamFlowcellData> FormatSpecificData;

    Layout(const boost::filesystem::path &baseCallsDirectory,
           const Format format,
           const FormatSpecificData &formatSpecificData,
           const unsigned laneNumberMax,
           const std::vector<unsigned> &barcodeCycles,
           const flowcell::ReadMetadataList &readMetadataList,
           const alignment::SeedMetadataList &seedMetadataList,
           const std::string &flowcellId);

    const boost::filesystem::path &getBaseCallsPath() const {return baseCallsPath_;}
    Format getFormat() const {return format_;}
    unsigned getLaneNumberMax() const {return laneNumberMax_;}
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
    const std::vector<unsigned> &getDataCycles() const {return dataCycles_;}

    unsigned getIndex() const {return index_;}
    void setIndex(unsigned index) {index_ = index;}

    template <Format format, typename AttributeTag> const typename AttributeTag::value_type& getAttribute(
        typename AttributeTag::value_type &result) const
    {
        BOOST_STATIC_ASSERT(sizeof(AttributeTag) == 0);
    }

    template <Format format, typename AttributeTag> typename AttributeTag::value_type getAttribute() const
    {
        typename AttributeTag::value_type ret;
        return getAttribute<format, AttributeTag>(ret);
    }

    template <Format format, typename AttributeTag> typename AttributeTag::value_type getLongestAttribute() const
    {
        BOOST_STATIC_ASSERT(sizeof(AttributeTag) == 0);
    }

    template <Format format, typename AttributeTag> void getLaneAttribute(
        const unsigned lane,
        typename AttributeTag::value_type& result) const
    {
        BOOST_STATIC_ASSERT(sizeof(AttributeTag) == 0);
    }

    template <Format format, typename AttributeTag> typename AttributeTag::value_type getLaneReadAttribute(
        const unsigned lane, const unsigned read) const
    {
        typename AttributeTag::value_type ret;
        return getLaneReadAttribute<format, AttributeTag>(lane, read, ret);
    }

    template <Format format, typename AttributeTag> const typename AttributeTag::value_type& getLaneReadAttribute(
        const unsigned lane, const unsigned read, typename AttributeTag::value_type &result) const
    {
        BOOST_STATIC_ASSERT(sizeof(AttributeTag) == 0);
    }

    template <Format format, typename AttributeTag> void getLaneCycleAttribute(
        const unsigned lane,
        const unsigned cycle,
        typename AttributeTag::value_type& result) const
    {
        BOOST_STATIC_ASSERT(sizeof(AttributeTag) == 0);
    }

    template <Format format, typename AttributeTag> void getLaneTileCycleAttribute(
        const unsigned lane,
        const unsigned tile,
        const unsigned cycle,
        typename AttributeTag::value_type& result) const
    {
        BOOST_STATIC_ASSERT(sizeof(AttributeTag) == 0);
    }

    template <Format format, typename AttributeTag> void getLaneTileAttribute(
        const unsigned lane,
        const unsigned tile,
        typename AttributeTag::value_type& result) const
    {
        BOOST_STATIC_ASSERT(sizeof(AttributeTag) == 0);
    }

private:
    boost::filesystem::path baseCallsPath_;
    Format format_;
    FormatSpecificData formatSpecificData_;
    unsigned laneNumberMax_;
    std::vector<unsigned> barcodeCycles_;
    std::string flowcellId_;
    // vector of sets at each lane number position.
    // Lowest lane number to date is 1, so position 0 must have an empty tile set
    typedef std::set<unsigned int> Tiles;
    typedef std::vector<Tiles >LaneTiles;
    LaneTiles laneTiles_;
    flowcell::ReadMetadataList readMetadataList_;
    alignment::SeedMetadataList seedMetadataList_;
    std::vector<unsigned> dataCycles_;
    unsigned index_;

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
        boost::bind(&ReadMetadataList::size, boost::bind(&Layout::getReadMetadataList, _1))<
        boost::bind(&ReadMetadataList::size, boost::bind(&Layout::getReadMetadataList, _2)))->getReadMetadataList().size();
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

inline unsigned getMaxCycleNumber(const FlowcellLayoutList &flowcellLayoutList)
{
    unsigned ret = 0;
    BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
    {
        ret = std::max(ret, getMaxCycleNumber(flowcell.getReadMetadataList()));
        if (!flowcell.getBarcodeCycles().empty())
        {
            ret = std::max(ret, *std::max_element(flowcell.getBarcodeCycles().begin(), flowcell.getBarcodeCycles().end()));
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

inline unsigned getMaxLaneNumber(const FlowcellLayoutList &flowcellLayoutList)
{
    unsigned ret = 0;
    BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
    {
        ret = std::max(ret, flowcell.getLaneNumberMax());
    }
    return ret;
}

template<Layout::Format format, typename AttributeTag>
typename AttributeTag::value_type getLongestAttribute(
    const std::vector<Layout> &flowcellLayoutList)
{
    typename AttributeTag::value_type ret;
    BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
    {
        if (format == flowcell.getFormat())
        {
            const typename AttributeTag::value_type attr = flowcell.getLongestAttribute<format, AttributeTag>();
            if (std::distance(ret.begin(), ret.end())< std::distance(attr.begin(), attr.end()))
            {
                ret = attr;
            }
        }
    }
    return ret;
}

template<Layout::Format format, typename AttributeTag>
typename AttributeTag::value_type getMaxAttribute(
    const std::vector<Layout> &flowcellLayoutList)
{
    typename AttributeTag::value_type ret(0);
    BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
    {
        if (format == flowcell.getFormat())
        {
            const typename AttributeTag::value_type attr = flowcell.getAttribute<format, AttributeTag>();
            ret = std::max(ret, attr);
        }
    }
    return ret;
}


inline std::ostream & operator << (std::ostream &os, const Layout &layout)
{
    os << "Layout("
              << layout.getFlowcellId() << ", "
              <<  "[";
    BOOST_FOREACH(const ReadMetadata &readMetadata, layout.getReadMetadataList())
    {
        os << "," << readMetadata;
    }
    return os
              << "])";
}

} // namespace flowcell
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_LAYOUT_HH
