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

std::pair<std::string, std::string> getSoftwareVersionFromConfig(const boost::filesystem::path &baseCallsPath)
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

Layout::Layout(const boost::filesystem::path &baseCallsDirectory,
        const Format format,
        const std::vector<unsigned> &barcodeCycles,
        const flowcell::ReadMetadataList &readMetadataList,
        const alignment::SeedMetadataList &seedMetadataList,
        const std::string &flowcellId)
     : baseCallsDirectory_(baseCallsDirectory)
     , format_(format)
     , barcodeCycles_(barcodeCycles)
     , flowcellId_(flowcellId)
     , readMetadataList_(readMetadataList)
     , seedMetadataList_(seedMetadataList)
     , allCycleNumbers_(flowcell::getAllCycleNumbers(readMetadataList_))
     , index_(0)
{
     if (Bcl == format_ || BclGz == format_)
     {
         softwareVersion_ = getSoftwareVersionFromConfig(baseCallsDirectory_);
     }
}


// TODO: there must be a better way to extract the path separator...
static const boost::filesystem::path a("a");
static const boost::filesystem::path aSlashA(a/a);
static const char directorySeparatorChar(aSlashA.string().at(1));

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
        else if (!strncmp(softwareVersion_.second.c_str(), "1.9.", 4) ||
            !strncmp(softwareVersion_.second.c_str(), "1.10.", 5) ||
            !strncmp(softwareVersion_.second.c_str(), "1.11.", 5) ||
            !strncmp(softwareVersion_.second.c_str(), "1.12.", 5) ||
            !strncmp(softwareVersion_.second.c_str(), "1.13.", 5) ||
            !strncmp(softwareVersion_.second.c_str(), "1.14.", 5) ||
            !strncmp(softwareVersion_.second.c_str(), "1.15.", 5) ||
            !strncmp(softwareVersion_.second.c_str(), "1.16.", 5) ||
            !strncmp(softwareVersion_.second.c_str(), "1.17.", 5))
        {
            pathInternalStringRef.append(laneFolder).append(filterFileName);
        }
        else
        {
            pathInternalStringRef.append(filterFileName);
        }
    }
    else
    {
        pathInternalStringRef.append(laneFolder).append(filterFileName);
    }
}


void Layout::getBclFilePath(
    const unsigned tile,
    const unsigned lane,
    const boost::filesystem::path &baseCallsPath,
    const unsigned cycle,
    const TileMetadata::Compression compression,
    boost::filesystem::path &result)
{
    ISAAC_ASSERT_MSG(lane <= maxLaneNumber_, "Lane number too big: " << lane);
    ISAAC_ASSERT_MSG(tile <= maxTileNumber_, "Tile number must not exceeed 4 digits");

    ISAAC_ASSERT_MSG(cycle <= maxCycleNumber_, "Cycle number should not exceeed 4 digits");
    // Warning: all this mad code below is to avoid memory allocations during path formatting.
    // the result is expected to be pre-sized, else allocations will occur as usual.
    char laneFolder[100];
    const int laneFolderLength = snprintf(laneFolder, sizeof(laneFolder), "%cL%03d", directorySeparatorChar, lane);

    char cycleFolder[100];
    const int cycleFolderLength = snprintf(cycleFolder, sizeof(cycleFolder), "%cC%d.1", directorySeparatorChar, cycle);

    char bclFileName[100];
    const int bclFileNameLength = snprintf(bclFileName, sizeof(bclFileName),
                                          (TileMetadata::GzCompression == compression ? "%cs_%d_%d.bcl.gz" : "%cs_%d_%d.bcl"),
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
