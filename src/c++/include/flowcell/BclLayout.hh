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
 ** \file BclLayout.hh
 **
 ** Specialization of Layout for bcl or bcl-gz flowcell.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_FLOWCELL_BCL_GZ_LAYOUT_HH
#define iSAAC_FLOWCELL_BCL_GZ_LAYOUT_HH

#include "flowcell/BclLayoutAttributes.hh"
#include "flowcell/Layout.hh"

namespace isaac
{
namespace flowcell
{

namespace bcl
{
    static const unsigned LANE_NUMBER_MAX = 8;
    static const unsigned TILE_NUMBER_MAX = 99999;
    static const unsigned CYCLE_NUMBER_MAX = 9999;
};

template<>
void Layout::getLaneTileCycleAttribute<Layout::Bcl, BclFilePathAttributeTag>(
    const unsigned lane, const unsigned tile, const unsigned cycle, boost::filesystem::path &result) const;

template<>
void Layout::getLaneTileAttribute<Layout::Bcl, FiltersFilePathAttributeTag>(
    const unsigned lane, const unsigned tile, boost::filesystem::path &result) const;

template<>
void Layout::getLaneTileAttribute<Layout::Bcl, PositionsFilePathAttributeTag>(
    const unsigned lane, const unsigned tile, boost::filesystem::path &result) const;

template<>
inline boost::filesystem::path Layout::getLongestAttribute<Layout::Bcl, BclFilePathAttributeTag>() const
{
    boost::filesystem::path cycleBciFilePath;
    getLaneTileCycleAttribute<Layout::Bcl, BclFilePathAttributeTag>(
        bcl::LANE_NUMBER_MAX, bcl::TILE_NUMBER_MAX, bcl::CYCLE_NUMBER_MAX, cycleBciFilePath);
    return cycleBciFilePath;
}

template<>
inline boost::filesystem::path Layout::getLongestAttribute<Layout::Bcl, FiltersFilePathAttributeTag>() const
{
    boost::filesystem::path filtersFilePath;
    getLaneTileAttribute<Layout::Bcl, FiltersFilePathAttributeTag>(bcl::LANE_NUMBER_MAX, bcl::TILE_NUMBER_MAX, filtersFilePath);
    return filtersFilePath;
}

template<>
boost::filesystem::path Layout::getLongestAttribute<Layout::Bcl, PositionsFilePathAttributeTag>() const;

} // namespace flowcell
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_BCL_GZ_LAYOUT_HH
