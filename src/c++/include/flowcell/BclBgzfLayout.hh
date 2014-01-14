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
 ** \file BclBgzfLayout.hh
 **
 ** Specialization of Layout for bcl-bgzf flowcell.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_FLOWCELL_BCL_BGZF_LAYOUT_HH
#define iSAAC_FLOWCELL_BCL_BGZF_LAYOUT_HH

#include "flowcell/BclLayoutAttributes.hh"
#include "flowcell/Layout.hh"

namespace isaac
{
namespace flowcell
{

namespace bclBgzf
{
    static const unsigned LANE_NUMBER_MAX = 8;
    static const unsigned CYCLE_NUMBER_MAX = 9999;
};

struct BciFilePathAttributeTag
{
    typedef boost::filesystem::path value_type;
    friend std::ostream &operator << (std::ostream &os, const BciFilePathAttributeTag &tag){return os << "BciFilePathAttributeTag";}
};

struct TilesPerLaneMaxAttributeTag
{
    typedef unsigned value_type;
    friend std::ostream &operator << (std::ostream &os, const TilesPerLaneMaxAttributeTag &tag){return os << "TilesPerLaneMaxAttributeTag";}
};

template<>
void Layout::getLaneCycleAttribute<Layout::BclBgzf, BclFilePathAttributeTag>(
    const unsigned lane, const unsigned cycle, boost::filesystem::path &result) const;

template<>
void Layout::getLaneAttribute<Layout::BclBgzf, FiltersFilePathAttributeTag>(
    const unsigned lane, boost::filesystem::path &result) const;

template<>
void Layout::getLaneAttribute<Layout::BclBgzf, PositionsFilePathAttributeTag>(
    const unsigned lane, boost::filesystem::path &result) const;

template<>
void Layout::getLaneCycleAttribute<Layout::BclBgzf, BciFilePathAttributeTag>(
    const unsigned lane, const unsigned cycle, boost::filesystem::path &result) const;

template<>
void Layout::getLaneAttribute<Layout::BclBgzf, BciFilePathAttributeTag>(
    const unsigned lane, boost::filesystem::path &result) const;

template<>
const unsigned& Layout::getAttribute<Layout::BclBgzf, TilesPerLaneMaxAttributeTag>(
    unsigned &result) const;

template<>
inline boost::filesystem::path Layout::getLongestAttribute<Layout::BclBgzf, BciFilePathAttributeTag>() const
{
    boost::filesystem::path cycleBciFilePath;
    getLaneCycleAttribute<Layout::BclBgzf, BciFilePathAttributeTag>(bclBgzf::LANE_NUMBER_MAX, bclBgzf::CYCLE_NUMBER_MAX, cycleBciFilePath);
    return cycleBciFilePath;
}

template<>
inline boost::filesystem::path Layout::getLongestAttribute<Layout::BclBgzf, BclFilePathAttributeTag>() const
{
    boost::filesystem::path cycleBclFilePath;
    getLaneCycleAttribute<Layout::BclBgzf, BclFilePathAttributeTag>(bclBgzf::LANE_NUMBER_MAX, bclBgzf::CYCLE_NUMBER_MAX, cycleBclFilePath);
    return cycleBclFilePath;
}

template<>
inline boost::filesystem::path Layout::getLongestAttribute<Layout::BclBgzf, FiltersFilePathAttributeTag>() const
{
    boost::filesystem::path filtersFilePath;
    getLaneAttribute<Layout::BclBgzf, FiltersFilePathAttributeTag>(bclBgzf::LANE_NUMBER_MAX, filtersFilePath);
    return filtersFilePath;
}

template<>
inline boost::filesystem::path Layout::getLongestAttribute<Layout::BclBgzf, PositionsFilePathAttributeTag>() const
{
    boost::filesystem::path positionsFilePath;
    getLaneAttribute<Layout::BclBgzf, PositionsFilePathAttributeTag>(bclBgzf::LANE_NUMBER_MAX, positionsFilePath);
    return positionsFilePath;
}

} // namespace flowcell
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_BCL_BGZF_LAYOUT_HH
