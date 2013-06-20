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

#include "basecalls/ConfigXml.hh"
#include "flowcell/Layout.hh"

namespace isaac
{
namespace flowcell
{

static std::pair<std::string, std::string> getSoftwareVersionFromConfig(const boost::filesystem::path &baseCallsPath)
{
    const boost::filesystem::path basecallsConfigXml(baseCallsPath / "config.xml");
    std::ifstream is(basecallsConfigXml.string().c_str());
    if (!is) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open basecalls config file " + basecallsConfigXml.string()));
    }
    basecalls::ConfigXml cfg;
    if (!(is >> cfg)) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to read from basecalls config file " + basecallsConfigXml.string()));
    }
    return cfg.getSoftwareVersion();
}

static std::pair<unsigned, unsigned> parseMajorMinor(const std::string &version)
{
    std::pair<unsigned, unsigned> ret;
    const std::size_t firstDotPos = version.find('.');
    if (std::string::npos == firstDotPos)
    {
        BOOST_THROW_EXCEPTION(common::UnsupportedVersionException(
            "Incorrect RTA software version format. Expected <major>.<minor>, got: " + version));
    }

    const std::string major = version.substr(0, firstDotPos);
    ret.first = boost::lexical_cast<unsigned>(major);
    const std::size_t secondDotPos = version.find('.', firstDotPos + 1);
    const std::string minor = (std::string::npos == secondDotPos) ?
        version.substr(firstDotPos + 1) : version.substr(firstDotPos + 1, secondDotPos - firstDotPos - 1);

    ret.second = boost::lexical_cast<unsigned>(minor);

    return ret;
}


Layout::Layout(const boost::filesystem::path &baseCallsDirectory,
        const Format format,
        const std::vector<unsigned> &barcodeCycles,
        const flowcell::ReadMetadataList &readMetadataList,
        const alignment::SeedMetadataList &seedMetadataList,
        const std::string &flowcellId)
     : baseCallsPath_(baseCallsDirectory)
     , format_(format)
     , barcodeCycles_(barcodeCycles)
     , flowcellId_(flowcellId)
     , readMetadataList_(readMetadataList)
     , seedMetadataList_(seedMetadataList)
     , allCycleNumbers_(flowcell::getAllCycleNumbers(readMetadataList_))
     , index_(0)
     , softwareMajorMinor_(std::make_pair(0U,0U))
{
     if (Bcl == format_ || BclGz == format_)
     {
         softwareVersion_ = getSoftwareVersionFromConfig(baseCallsPath_);
         softwareMajorMinor_ = parseMajorMinor(softwareVersion_.second);
     }
}


// TODO: there must be a better way to extract the path separator...
static const boost::filesystem::path a("a");
static const boost::filesystem::path aSlashA(a/a);
static const char directorySeparatorChar(aSlashA.string().at(1));

void Layout::getBamFilePath(
    const boost::filesystem::path &baseCallsPath,
    boost::filesystem::path &result)
{
    // boost 1.46 implementation of filesystem::path is coded to instantiate an std::string
    // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());
    pathInternalStringRef = baseCallsPath.string().c_str();
}

std::size_t Layout::getBamFileSize() const
{
    boost::filesystem::path filePath;
    getBamFilePath(filePath);
    return boost::filesystem::file_size(filePath);
}

void Layout::getFastqFilePath(
    const unsigned read,
    const unsigned lane,
    const boost::filesystem::path &baseCallsPath,
    const bool compressed,
    boost::filesystem::path &result)
{
    ISAAC_ASSERT_MSG(lane <= maxLaneNumber_, "Lane number " << lane << " must not exceed " << maxLaneNumber_);
    ISAAC_ASSERT_MSG(read <= maxReadNumber_, "Read number " << read << " must not exceed " << maxReadNumber_);

    // Warning: all this mad code below is to avoid memory allocations during path formatting.
    // the result is expected to be pre-sized, else allocations will occur as usual.
    char fastqFileName[100];
    const int fastqFileNameLength = snprintf(fastqFileName, sizeof(fastqFileName),
                                          (compressed ?
                                              "%clane%d_read%d.fastq.gz" : "%clane%d_read%d.fastq"),
                                          directorySeparatorChar, lane, read);

    // boost 1.46 implementation of filesystem::path is coded to instantiate an std::string
    // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());
    pathInternalStringRef = baseCallsPath.string().c_str();
    pathInternalStringRef.append(fastqFileName, fastqFileName + fastqFileNameLength);
}

void Layout::getFiltersFilePath(
    const unsigned tile,
    const unsigned lane,
    const boost::filesystem::path &baseCallsPath,
    boost::filesystem::path &result) const
{
    ISAAC_ASSERT_MSG(lane <= maxLaneNumber_, "Lane number too big: " << lane);
    ISAAC_ASSERT_MSG(tile <= maxTileNumber_, "Tile number must not exceed 4 digits : " << tile);

    // Warning: all this mad code below is to avoid memory allocations during path formatting.
    // the result is expected to be pre-sized, else allocations will occur as usual.
    char laneFolder[100];
    sprintf(laneFolder, "%cL%03d", directorySeparatorChar, lane);

    // boost 1.46 implementation of filesystem::path is coded to instantiated std::string
    // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());

    char filterFileName[100];
    sprintf(filterFileName, "%cs_%d_%04d.filter", directorySeparatorChar, lane, tile);

    pathInternalStringRef = baseCallsPath.c_str();
    if (softwareVersion_.first == "RTA")
    {
        if (!strncmp(softwareVersion_.second.c_str(), "1.8.", 4))
        {
            pathInternalStringRef.append(filterFileName);
        }
        else
        {
            if (softwareMajorMinor_.first > 1 || (softwareMajorMinor_.first == 1 && softwareMajorMinor_.second >= 9))
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

void Layout::getPositionsFilePath(
    const unsigned tile,
    const unsigned lane,
    const boost::filesystem::path &baseCallsPath,
    boost::filesystem::path &result) const
{
    ISAAC_ASSERT_MSG(lane <= maxLaneNumber_, "Lane number too big: " << lane);
    ISAAC_ASSERT_MSG(tile <= maxTileNumber_, "Tile number must not exceed 4 digits : " << tile);

    // Warning: all this mad code below is to avoid memory allocations during path formatting.
    // the result is expected to be pre-sized, else allocations will occur as usual.
    char laneFolder[100];
    // assuming Intensities folder is one level anove BaseCalls folder
    sprintf(laneFolder, "%c..%cL%03d", directorySeparatorChar, directorySeparatorChar, lane);

    // boost 1.46 implementation of filesystem::path is coded to instantiated std::string
    // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());

    char filterFileName[100];
    sprintf(filterFileName, "%cs_%d_%04d.clocs", directorySeparatorChar, lane, tile);

    pathInternalStringRef = baseCallsPath.c_str();

    pathInternalStringRef.append(laneFolder).append(filterFileName);
}


void Layout::getBclFilePath(
    const unsigned tile,
    const unsigned lane,
    const boost::filesystem::path &baseCallsPath,
    const unsigned cycle,
    const bool compressed,
    boost::filesystem::path &result)
{
    ISAAC_ASSERT_MSG(lane <= maxLaneNumber_, "Lane number too big: " << lane);
    ISAAC_ASSERT_MSG(tile <= maxTileNumber_, "Tile number must not exceeed 4 digits. got: " << tile);

    ISAAC_ASSERT_MSG(cycle <= maxCycleNumber_, "Cycle number should not exceeed 4 digits");
    // Warning: all this mad code below is to avoid memory allocations during path formatting.
    // the result is expected to be pre-sized, else allocations will occur as usual.
    char laneFolder[100];
    const int laneFolderLength = snprintf(laneFolder, sizeof(laneFolder), "%cL%03d", directorySeparatorChar, lane);

    char cycleFolder[100];
    const int cycleFolderLength = snprintf(cycleFolder, sizeof(cycleFolder), "%cC%d.1", directorySeparatorChar, cycle);

    char bclFileName[100];
    const int bclFileNameLength = snprintf(bclFileName, sizeof(bclFileName),
                                          (compressed ? "%cs_%d_%d.bcl.gz" : "%cs_%d_%d.bcl"),
                                          directorySeparatorChar, lane, tile);

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


} // namespace flowcell
} // namespace isaac
