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
 ** \file Layout.cpp
 **
 ** Flowcell file locations and such
 **
 ** \author Roman Petrovski
 **/

#include "common/FileSystem.hh"
#include "flowcell/BclLayout.hh"

namespace isaac
{
namespace flowcell
{

static void getPositionsFilePath(
    const boost::filesystem::path &baseCallsPath,
    const unsigned lane,
    const unsigned tile,
    const bool patternedFlowcell,
    const bool forceClocs,
    boost::filesystem::path &result)
{
    ISAAC_ASSERT_MSG(lane <= bcl::LANE_NUMBER_MAX, "Lane number " << lane << " must not exceed " << bcl::LANE_NUMBER_MAX);
    ISAAC_ASSERT_MSG(tile <= bcl::TILE_NUMBER_MAX, "Tile number must not exceed 5 digits : " << tile);

    if (patternedFlowcell)
    {
        result = baseCallsPath.c_str();
        result /= "..";
        result /= forceClocs ? "s.clocs" : "s.locs";
    }
    else
    {
        // Warning: all this mad code below is to avoid memory allocations during path formatting.
        // the result is expected to be pre-sized, else allocations will occur as usual.
        char laneFolder[100];
        // assuming Intensities folder is one level anove BaseCalls folder
        sprintf(laneFolder, "%c..%cL%03d", common::getDirectorySeparatorChar(), common::getDirectorySeparatorChar(), lane);

        // boost 1.46 implementation of filesystem::path is coded to instantiated std::string
        // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
    //    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());

        char filterFileName[100];
        // most instruments are configured to produce clocs, used it by default
        sprintf(filterFileName, "%cs_%d_%04d.clocs", common::getDirectorySeparatorChar(), lane, tile);

        result = baseCallsPath.c_str();
        result /= laneFolder;
        result /= filterFileName;

        if (!forceClocs && !boost::filesystem::exists(result))
        {
            // assume it's uncompressed positions then...
            sprintf(filterFileName, "%cs_%d_%04d.locs", common::getDirectorySeparatorChar(), lane, tile);

            result = baseCallsPath.c_str();
            result /= laneFolder;
            result /= filterFileName;
        }
    }
}

static void getBclFilePath(
    const boost::filesystem::path &baseCallsPath,
    const unsigned lane,
    const unsigned tile,
    const unsigned cycle,
    const bool compressed,
    boost::filesystem::path &result)
{
    ISAAC_ASSERT_MSG(lane <= bcl::LANE_NUMBER_MAX, "Lane number " << lane << " must not exceed " << bcl::LANE_NUMBER_MAX);
    ISAAC_ASSERT_MSG(tile <= bcl::TILE_NUMBER_MAX, "Tile number must not exceeed 5 digits. got: " << tile);

    ISAAC_ASSERT_MSG(cycle <= bcl::CYCLE_NUMBER_MAX, "Cycle number should not exceeed 4 digits");
    // Warning: all this mad code below is to avoid memory allocations during path formatting.
    // the result is expected to be pre-sized, else allocations will occur as usual.
    char laneFolder[100];
    const int laneFolderLength = snprintf(laneFolder, sizeof(laneFolder), "%cL%03d", common::getDirectorySeparatorChar(), lane);

    char cycleFolder[100];
    const int cycleFolderLength = snprintf(cycleFolder, sizeof(cycleFolder), "%cC%d.1", common::getDirectorySeparatorChar(), cycle);

    char bclFileName[100];
    const int bclFileNameLength = snprintf(bclFileName, sizeof(bclFileName),
                                          (compressed ? "%cs_%d_%d.bcl.gz" : "%cs_%d_%d.bcl"),
                                          common::getDirectorySeparatorChar(), lane, tile);

    // boost 1.46 implementation of filesystem::path is coded to instantiate an std::string
    // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());
    pathInternalStringRef = baseCallsPath.string().c_str();
    pathInternalStringRef.append(laneFolder, laneFolder + laneFolderLength);
    pathInternalStringRef.append(cycleFolder, cycleFolder + cycleFolderLength);
    pathInternalStringRef.append(bclFileName, bclFileName + bclFileNameLength);

//    std::cerr << "formatted " << result << " out of " << laneFolder << ","
//        << cycleFolder << "," << bclFileName << "\n";
}



template<>
void Layout::getLaneTileCycleAttribute<Layout::Bcl, BclFilePathAttributeTag>(
    const unsigned lane, const unsigned tile, const unsigned cycle, boost::filesystem::path &result) const
{
    ISAAC_ASSERT_MSG(Bcl == format_, BclFilePathAttributeTag() << " is only allowed for bcl and bcl-gz flowcells");

    const BclFlowcellData &data = boost::get<BclFlowcellData>(formatSpecificData_);
    return flowcell::getBclFilePath(getBaseCallsPath(), lane, tile, cycle, data.compressed_, result);
}

template<>
void Layout::getLaneTileAttribute<Layout::Bcl, FiltersFilePathAttributeTag>(
    const unsigned lane, const unsigned tile, boost::filesystem::path &result) const
{
    ISAAC_ASSERT_MSG(Bcl == format_, FiltersFilePathAttributeTag() << " is only allowed for bcl and bcl-gz flowcells");

    ISAAC_ASSERT_MSG(lane <= bcl::LANE_NUMBER_MAX, "Lane number " << lane << " must not exceed " << bcl::LANE_NUMBER_MAX);
    ISAAC_ASSERT_MSG(tile <= bcl::TILE_NUMBER_MAX, "Tile number must not exceed 5 digits : " << tile);

    // Warning: all this mad code below is to avoid memory allocations during path formatting.
    // the result is expected to be pre-sized, else allocations will occur as usual.
    char laneFolder[100];
    sprintf(laneFolder, "%cL%03d", common::getDirectorySeparatorChar(), lane);

    // boost 1.46 implementation of filesystem::path is coded to instantiated std::string
    // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());

    char filterFileName[100];
    sprintf(filterFileName, "%cs_%d_%04d.filter", common::getDirectorySeparatorChar(), lane, tile);

    pathInternalStringRef = getBaseCallsPath().c_str();
    const BclFlowcellData &data = boost::get<BclFlowcellData>(formatSpecificData_);

    if (data.softwareVersion_.first == "RTA")
    {
        if (!strncmp(data.softwareVersion_.second.c_str(), "1.8.", 4))
        {
            pathInternalStringRef.append(filterFileName);
        }
        else
        {
            if (data.softwareMajorMinor_.first > 1 || (data.softwareMajorMinor_.first == 1 && data.softwareMajorMinor_.second >= 9))
            {
                pathInternalStringRef.append(laneFolder).append(filterFileName);
            }
            else
            {
                pathInternalStringRef.append(filterFileName);
            }
        }
    }
    else
    {
        pathInternalStringRef.append(laneFolder).append(filterFileName);
    }
}

template<>
void Layout::getLaneTileAttribute<Layout::Bcl, PositionsFilePathAttributeTag>(
    const unsigned lane, const unsigned tile, boost::filesystem::path &result) const
{
    ISAAC_ASSERT_MSG(Bcl == format_, PositionsFilePathAttributeTag() << " is only allowed for bcl and bcl-gz flowcells");

    const BclFlowcellData &data = boost::get<BclFlowcellData>(formatSpecificData_);

    return flowcell::getPositionsFilePath(getBaseCallsPath(), lane, tile, data.patternedFlowcell_, false, result);
}

template<>
boost::filesystem::path Layout::getLongestAttribute<Layout::Bcl, PositionsFilePathAttributeTag>() const
{
    const BclFlowcellData &data = boost::get<BclFlowcellData>(formatSpecificData_);
    boost::filesystem::path positionsFilePath;
    flowcell::getPositionsFilePath(
        getBaseCallsPath(), bcl::LANE_NUMBER_MAX, bcl::TILE_NUMBER_MAX,
        data.patternedFlowcell_,
        // .clocs is longer than .locs
        true,
        positionsFilePath);
    return positionsFilePath;
}

} // namespace flowcell
} // namespace isaac
