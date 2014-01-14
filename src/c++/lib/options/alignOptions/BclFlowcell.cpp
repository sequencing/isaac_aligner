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
 ** \file BclFlowcell.hh
 **
 ** Generate flowcell object out of BaseCalls/config.xml
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_attr_cast.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix1_statements.hpp>
#include <boost/spirit/home/phoenix/statement/if.hpp>
#include <boost/spirit/home/phoenix/statement/throw.hpp>
#include <boost/spirit/home/phoenix/object/static_cast.hpp>
#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/spirit/home/phoenix/bind.hpp>
//#include <boost/spirit/home/phoenix/object/construct.hpp>
//#include <boost/spirit/home/support/unused.hpp>

#include <boost/algorithm/string/regex.hpp>

#include "rta/ConfigXml.hh"
#include "rta/RunInfoXml.hh"

#include "alignOptions/UseBasesMaskOption.hh"

#include "BclFlowcell.hh"
#include "SeedDescriptorOption.hh"


namespace isaac
{
namespace options
{
namespace alignOptions
{

struct ParsedBaseCallsConfig
{
    std::string flowcellId_;
    struct ReadInfo
    {
        ReadInfo(const unsigned firstCycle, const unsigned lastCycle) :
            firstCycle_(firstCycle), lastCycle_(lastCycle){}
        unsigned firstCycle_;
        unsigned lastCycle_;
    };
    std::vector<ReadInfo> readInfo_;
    std::vector<unsigned> lanes_;
    std::vector<std::vector<unsigned> > laneTiles_;
    boost::filesystem::path baseCallsPath_;

    flowcell::BclFlowcellData bclFlowcellData_;
};


static std::pair<std::string, std::string> getSoftwareVersionFromConfig(const boost::filesystem::path &baseCallsPath)
{
    const boost::filesystem::path basecallsConfigXml(baseCallsPath / "config.xml");
    std::ifstream is(basecallsConfigXml.string().c_str());
    if (!is) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open basecalls config file " + basecallsConfigXml.string()));
    }
    rta::ConfigXml cfg;
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

static ParsedBaseCallsConfig parseBasecallsConfigXml(
    const bool compressed,
    const boost::filesystem::path &basecallsConfigXml)
{
    std::ifstream is(basecallsConfigXml.string().c_str());
    if (!is) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open basecalls config file " + basecallsConfigXml.string()));
    }
    rta::ConfigXml cfg;
    if (!(is >> cfg)) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to read from basecalls config file " + basecallsConfigXml.string()));
    }

    ParsedBaseCallsConfig ret;
    ret.flowcellId_ = cfg.getFlowcellId();
    BOOST_FOREACH(const rta::ConfigXml::RunParametersRead &runParametersRead, cfg.getRunParametersReads())
    {
        ret.readInfo_.push_back(ParsedBaseCallsConfig::ReadInfo(runParametersRead.firstCycle_, runParametersRead.lastCycle_));
    }
    ret.lanes_ = cfg.getLanes();
    BOOST_FOREACH(const unsigned lane, ret.lanes_)
    {
        ret.laneTiles_.resize(std::max<std::size_t>(ret.laneTiles_.size(), lane + 1));
        ret.laneTiles_.at(lane) = cfg.getTiles(lane);
    }

    ret.baseCallsPath_ = basecallsConfigXml.parent_path();

    ret.bclFlowcellData_.softwareVersion_ = getSoftwareVersionFromConfig(ret.baseCallsPath_);
    ret.bclFlowcellData_.softwareMajorMinor_ = parseMajorMinor(ret.bclFlowcellData_.softwareVersion_.second);
    ret.bclFlowcellData_.compressed_ = compressed;
    ret.bclFlowcellData_.patternedFlowcell_ = boost::filesystem::exists(ret.baseCallsPath_.parent_path() / "s.locs");
    ret.bclFlowcellData_.tilesPerLaneMax_ = std::max<unsigned>(ret.bclFlowcellData_.tilesPerLaneMax_, ret.laneTiles_.size());

    return ret;
}

static ParsedBaseCallsConfig parseRunInfoXml(
    const bool compressed,
    const boost::filesystem::path &basecallsConfigXml)
{
    std::ifstream is(basecallsConfigXml.string().c_str());
    if (!is) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open basecalls config file " + basecallsConfigXml.string()));
    }
    rta::RunInfoXml cfg;
    if (!(is >> cfg)) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to read from basecalls config file " + basecallsConfigXml.string()));
    }

    ParsedBaseCallsConfig ret;
    ret.flowcellId_ = cfg.getFlowcellId();
    unsigned currentCycle = 1;
    BOOST_FOREACH(const rta::RunInfoXml::ReadInfo &readInfo, cfg.getReadInfos())
    {
        ret.readInfo_.push_back(ParsedBaseCallsConfig::ReadInfo(currentCycle, currentCycle + readInfo.numberOfCycles_ - 1));
        currentCycle += readInfo.numberOfCycles_ ;
    }
    ret.lanes_ = cfg.getLanes();
    BOOST_FOREACH(const unsigned lane, ret.lanes_)
    {
        ret.laneTiles_.resize(std::max<std::size_t>(ret.laneTiles_.size(), lane + 1));
        ret.laneTiles_.at(lane) = cfg.getTiles(lane);
    }

    ret.baseCallsPath_ = basecallsConfigXml.parent_path() / "Data" / "Intensities" / "BaseCalls";
    ret.bclFlowcellData_.compressed_ = compressed;
    ret.bclFlowcellData_.patternedFlowcell_ = boost::filesystem::exists(ret.baseCallsPath_.parent_path() / "s.locs");
    ret.bclFlowcellData_.tilesPerLaneMax_ = std::max<unsigned>(ret.bclFlowcellData_.tilesPerLaneMax_, ret.laneTiles_.size());
    return ret;
}

static ParsedBaseCallsConfig parseBaseCallsMetadata(
    const flowcell::Layout::Format format,
    const bool compressed,
    boost::filesystem::path baseCallsPath)
{

    if (boost::filesystem::is_directory(baseCallsPath))
    {
        baseCallsPath /= "config.xml";
        if (!boost::filesystem::exists(baseCallsPath))
        {
            const boost::format message =
                boost::format("\n   *** File not found: %s. config.xml must exist if --base-calls points to a folder. "
                    "Otherwise please supply path to RunInfo.xml ***\n") %
                baseCallsPath.string();
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }
        // Assuming this is BaseCalls folder path, try opening config.xml for backward compatibility
        return parseBasecallsConfigXml(compressed, baseCallsPath);
    }
    return parseRunInfoXml(compressed, baseCallsPath);
}

flowcell::Layout BclFlowcell::createFilteredFlowcell(
    const bool detectSimpleIndels,
    const std::string &tilesFilter,
    const boost::filesystem::path &baseCallsPath,
    const flowcell::Layout::Format format,
    const bool compressed,
    const unsigned laneNumberMax,
    std::string useBasesMask,
    const std::string &seedDescriptor,
    const unsigned seedLength,
    const reference::ReferenceMetadataList &referenceMetadataList,
    unsigned &firstPassSeeds)
{
    const ParsedBaseCallsConfig cfg = parseBaseCallsMetadata(format, compressed, baseCallsPath);

    using boost::phoenix::bind;
    using boost::phoenix::arg_names::_1;

    std::vector<unsigned int> readLengths;
    std::transform(cfg.readInfo_.begin(), cfg.readInfo_.end(),
                   std::back_inserter(readLengths),
                   bind(&ParsedBaseCallsConfig::ReadInfo::lastCycle_, _1) -
                       bind(&ParsedBaseCallsConfig::ReadInfo::firstCycle_, _1) + 1);

    // TODO: this is guessing for the poor. Implement the proper one based on RunInfo.xml. config.xml does not contain
    // proper information about second barcode read in RTA 1.13.46.0
    if ("default" == useBasesMask)
    {
        if (readLengths.size() == 1)
        {
            useBasesMask = "y*n";
        }
        else if (readLengths.size() == 2)
        {
            useBasesMask = "y*n,y*n";
        }
        else if (readLengths.size() == 3)
        {
            useBasesMask = "y*n,i*n,y*n";
        }
        else if (readLengths.size() == 4)
        {
            // in zebra masking the last barcode cycles must not happen
            useBasesMask = "y*n,i*,i*,y*n";
        }
        else
        {
            const boost::format message =
                boost::format("\n   *** Could not guess the use-bases-mask for '%s', please supply the explicit value ***\n") %
                baseCallsPath.string();
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }

    }

    std::vector<unsigned int> readFirstCycles;
    std::transform(cfg.readInfo_.begin(), cfg.readInfo_.end(),
                   std::back_inserter(readFirstCycles),
                   bind(&ParsedBaseCallsConfig::ReadInfo::firstCycle_, _1));


    ParsedUseBasesMask parsedUseBasesMask = parseUseBasesMask(readFirstCycles, readLengths, seedLength, useBasesMask, baseCallsPath);
    std::vector<unsigned> barcodeCycles;
    BOOST_FOREACH(const flowcell::ReadMetadata &barcodeRead, parsedUseBasesMask.indexReads_)
    {
        barcodeCycles.insert(barcodeCycles.end(), barcodeRead.getCycles().begin(), barcodeRead.getCycles().end());
    }

    const alignment::SeedMetadataList seedMetadataList =
        parseSeedDescriptor(detectSimpleIndels, parsedUseBasesMask.dataReads_, seedDescriptor, seedLength, firstPassSeeds);

    flowcell::Layout fc(cfg.baseCallsPath_,
                        format,
                        cfg.bclFlowcellData_,
                        laneNumberMax,
                        barcodeCycles,
                        parsedUseBasesMask.dataReads_,
                        seedMetadataList, cfg.flowcellId_);

    std::string regexString(tilesFilter);
    std::replace(regexString.begin(), regexString.end(), ',', '|');
    boost::regex re(regexString);
    BOOST_FOREACH(const unsigned int lane, cfg.lanes_)
    {
        BOOST_FOREACH(const unsigned int tile, cfg.laneTiles_.at(lane))
        {
            std::string tileString((boost::format("s_%d_%04d") % lane % tile).str());
            if (boost::regex_search(tileString, re))
            {
                fc.addTile(lane, tile);
            }
        }
    }

    if (fc.getLaneIds().empty())
    {
        const boost::format message = boost::format("\n   *** Could not find any tiles matching the '%s' in: %s ***\n") %
            tilesFilter % baseCallsPath;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    return fc;
}


} // namespace alignOptions
} // namespace option
} // namespace isaac
