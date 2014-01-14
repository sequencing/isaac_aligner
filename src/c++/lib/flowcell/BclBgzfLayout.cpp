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
 ** \file Layout.cpp
 **
 ** Flowcell file locations and such
 **
 ** \author Roman Petrovski
 **/

#include "common/FileSystem.hh"
#include "flowcell/BclBgzfLayout.hh"

namespace isaac
{
namespace flowcell
{

static void getFiltersFilePath(
    const boost::filesystem::path &baseCallsPath,
    const unsigned lane,
    boost::filesystem::path &result)
{
    ISAAC_ASSERT_MSG(lane <= bclBgzf::LANE_NUMBER_MAX, "Lane number " << lane << " must not exceed " << bclBgzf::LANE_NUMBER_MAX);

    // Warning: all this mad code below is to avoid memory allocations during path formatting.
    // the result is expected to be pre-sized, else allocations will occur as usual.
    char laneFolder[100];
    sprintf(laneFolder, "%cL%03d", common::getDirectorySeparatorChar(), lane);

    // boost 1.46 implementation of filesystem::path is coded to instantiated std::string
    // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());

    char filterFileName[100];
    sprintf(filterFileName, "%cs_%d.filter", common::getDirectorySeparatorChar(), lane);

    pathInternalStringRef = baseCallsPath.c_str();
    pathInternalStringRef.append(laneFolder).append(filterFileName);
}

static void getPositionsFilePath(
    const boost::filesystem::path &baseCallsPath,
    const unsigned lane,
    boost::filesystem::path &result)
{
    ISAAC_ASSERT_MSG(lane <= bclBgzf::LANE_NUMBER_MAX, "Lane number " << lane << " must not exceed " << bclBgzf::LANE_NUMBER_MAX);

    // Warning: all this mad code below is to avoid memory allocations during path formatting.
    // the result is expected to be pre-sized, else allocations will occur as usual.
    char laneFolder[100];
    // assuming Intensities folder is one level anove BaseCalls folder
    sprintf(laneFolder, "%c..%cL%03d", common::getDirectorySeparatorChar(), common::getDirectorySeparatorChar(), lane);

    // boost 1.46 implementation of filesystem::path is coded to instantiated std::string
    // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
//    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());

    char filterFileName[100];
    sprintf(filterFileName, "%cs_%d.locs", common::getDirectorySeparatorChar(), lane);

    result = baseCallsPath.c_str();
    result /= laneFolder;
    result /= filterFileName;
}

static void getBclBgzfCycleFilePath(
    const boost::filesystem::path &baseCallsPath,
    const unsigned lane,
    const unsigned cycle,
    const bool bci,
    boost::filesystem::path &result)
{
    ISAAC_ASSERT_MSG(lane <= bclBgzf::LANE_NUMBER_MAX, "Lane number " << lane << " must not exceed " << bclBgzf::LANE_NUMBER_MAX);
    ISAAC_ASSERT_MSG(cycle <= bclBgzf::CYCLE_NUMBER_MAX, "Cycle number should not exceeed " <<
                     bclBgzf::CYCLE_NUMBER_MAX << " digits");
    // Warning: all this mad code below is to avoid memory allocations during path formatting.
    // the result is expected to be pre-sized, else allocations will occur as usual.
    char laneFolder[100];
    const int laneFolderLength = snprintf(laneFolder, sizeof(laneFolder), "%cL%03d", common::getDirectorySeparatorChar(), lane);

    char bclFileName[100];
    const int bclFileNameLength = snprintf(bclFileName, sizeof(bclFileName),
                                           bci ? "%c%04d.bcl.bgzf.bci":"%c%04d.bcl.bgzf", common::getDirectorySeparatorChar(), cycle);

    // boost 1.46 implementation of filesystem::path is coded to instantiate an std::string
    // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());
    pathInternalStringRef = baseCallsPath.c_str();
    pathInternalStringRef.append(laneFolder, laneFolder + laneFolderLength);
    pathInternalStringRef.append(bclFileName, bclFileName + bclFileNameLength);

//    std::cerr << "formatted " << result << " out of " << laneFolder << ","
//        << cycleFolder << "," << bclFileName << "\n";
}

static void getLaneBciFilePath(
    const boost::filesystem::path &baseCallsPath,
    const unsigned lane,
    boost::filesystem::path &result)
{
    ISAAC_ASSERT_MSG(lane <= bclBgzf::LANE_NUMBER_MAX, "Lane number " << lane << " must not exceed " << bclBgzf::LANE_NUMBER_MAX);

    // Warning: all this mad code below is to avoid memory allocations during path formatting.
    // the result is expected to be pre-sized, else allocations will occur as usual.
    char laneFolder[100];
    const int laneFolderLength = snprintf(laneFolder, sizeof(laneFolder), "%cL%03d", common::getDirectorySeparatorChar(), lane);

    char bciFileName[100];
    const int bciFileNameLength = snprintf(bciFileName, sizeof(bciFileName), "%cs_%d.bci", common::getDirectorySeparatorChar(), lane);

    // boost 1.46 implementation of filesystem::path is coded to instantiate an std::string
    // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());
    pathInternalStringRef = baseCallsPath.c_str();
    pathInternalStringRef.append(laneFolder, laneFolder + laneFolderLength);
    pathInternalStringRef.append(bciFileName, bciFileName + bciFileNameLength);
}

template<>
void Layout::getLaneCycleAttribute<Layout::BclBgzf, BclFilePathAttributeTag>(
    const unsigned lane, const unsigned cycle, boost::filesystem::path &result) const
{
    ISAAC_ASSERT_MSG(BclBgzf == format_, BclFilePathAttributeTag() << " is only allowed for bcl-bgzf flowcells");

    return getBclBgzfCycleFilePath(getBaseCallsPath(), lane, cycle, false, result);
}

template<>
void Layout::getLaneAttribute<Layout::BclBgzf, FiltersFilePathAttributeTag>(
    const unsigned lane, boost::filesystem::path &result) const
{
    ISAAC_ASSERT_MSG(BclBgzf == format_, FiltersFilePathAttributeTag() << " is only allowed for bcl-bgzf flowcells");

    return flowcell::getFiltersFilePath(getBaseCallsPath(), lane, result);
}

template<>
void Layout::getLaneAttribute<Layout::BclBgzf, PositionsFilePathAttributeTag>(
    const unsigned lane, boost::filesystem::path &result) const
{
    ISAAC_ASSERT_MSG(BclBgzf == format_, PositionsFilePathAttributeTag() << " is only allowed for bcl-bgzf flowcells");

    return flowcell::getPositionsFilePath(getBaseCallsPath(), lane, result);
}

template<>
void Layout::getLaneCycleAttribute<Layout::BclBgzf, BciFilePathAttributeTag>(
    const unsigned lane, const unsigned cycle, boost::filesystem::path &result) const
{
    ISAAC_ASSERT_MSG(BclBgzf == format_, BciFilePathAttributeTag() << " is only allowed for bcl-bgzf flowcells");

    return getBclBgzfCycleFilePath(getBaseCallsPath(), lane, cycle, true, result);
}

template<>
void Layout::getLaneAttribute<Layout::BclBgzf, BciFilePathAttributeTag>(
    const unsigned lane, boost::filesystem::path &result) const
{
    ISAAC_ASSERT_MSG(BclBgzf == format_, BciFilePathAttributeTag() << " is only allowed for bcl-bgzf flowcells");

    return getLaneBciFilePath(getBaseCallsPath(), lane, result);
}

template<>
const unsigned& Layout::getAttribute<Layout::BclBgzf, TilesPerLaneMaxAttributeTag>(
    unsigned &result) const
{
    ISAAC_ASSERT_MSG(BclBgzf == format_, TilesPerLaneMaxAttributeTag() << " is only allowed for bcl-bgzf flowcells");

    const BclFlowcellData &data = boost::get<BclFlowcellData>(formatSpecificData_);

    return data.tilesPerLaneMax_;
}

} // namespace flowcell
} // namespace isaac
