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

#include <boost/algorithm/string/regex.hpp>

#include "alignOptions/UseBasesMaskOption.hh"

#include "common/Threads.hpp"
#include "io/BamLoader.hh"

#include "BamFlowcell.hh"
#include "SeedDescriptorOption.hh"


namespace isaac
{
namespace options
{
namespace alignOptions
{

BamFlowcell::BamPath BamFlowcell::findBamPath(
    const boost::filesystem::path &baseCallsPath)
{
    BamPath ret;

    if (boost::filesystem::exists(baseCallsPath))
    {
        if (boost::filesystem::is_regular_file(baseCallsPath))
        {
            ret.path_ = baseCallsPath;
            ret.lane_ = 1;
        }
        else
        {
            BOOST_THROW_EXCEPTION(common::InvalidOptionException("Bam --base-calls must be a regular file. Got: "  + baseCallsPath.string()));
        }
    }
    else
    {
        BOOST_THROW_EXCEPTION(common::InvalidOptionException("Bam file does not exist: "  + baseCallsPath.string()));
    }

    return ret;
}

void nothing()
{
}

class MetadataParser
{
    BamFlowcellInfo &flowcellInfo_;
    bool pairednessKnown_;
    bool paired_;

public:
    MetadataParser(BamFlowcellInfo &flowcellInfo) : flowcellInfo_(flowcellInfo), pairednessKnown_(false), paired_(false)
    {
    }

    bool parseMetadata(
        const bool allowVariableReadLength,
        const bool allowMixedFlowcells,
        const bam::BamBlockHeader &block,
        const bool lastBlock)
    {
        const std::string flowcellId = parseFlowcellId(block);
//        ISAAC_THREAD_CERR << flowcellId << std::endl;
        if (flowcellInfo_.flowcellId_.empty())
        {
            flowcellInfo_.flowcellId_ = flowcellId;
        }
        else if (flowcellInfo_.flowcellId_ != flowcellId && !allowMixedFlowcells)
        {
            BOOST_THROW_EXCEPTION(common::InvalidOptionException("Multiple flowcells detected in the bam file: " +
                flowcellInfo_.flowcellId_ + " and " + flowcellId + ". Please specify an explicit --use-bases-mask to enable mixed flowcells."));
        }

        if (!pairednessKnown_)
        {
            paired_ = block.isPaired();
            pairednessKnown_ = true;
        }
        else if (paired_ != block.isPaired())
        {
            BOOST_THROW_EXCEPTION(common::InvalidOptionException("Mix of paired and single-ended data is not supported."));
        }

        ISAAC_ASSERT_MSG(pairednessKnown_, "It should be enough to see one segment to know if bam is paired or not");
        if (block.isReadOne())
        {
            if (!flowcellInfo_.readLengths_.first)
            {
                flowcellInfo_.readLengths_.first = block.getLSeq();
            }
            else if (!allowVariableReadLength && flowcellInfo_.readLengths_.first != unsigned(block.getLSeq()))
            {
                BOOST_THROW_EXCEPTION(common::InvalidOptionException("Mix of varying read lengths is not supported. Found: " +
                    boost::lexical_cast<std::string>(flowcellInfo_.readLengths_.first) + " and " + boost::lexical_cast<std::string>(block.getLSeq()) +
                    " for read 1"));
            }
        }
        else
        {
            if (!flowcellInfo_.readLengths_.second)
            {
                flowcellInfo_.readLengths_.second = block.getLSeq();
            }
            else if (!allowVariableReadLength && flowcellInfo_.readLengths_.second != unsigned(block.getLSeq()))
            {
                BOOST_THROW_EXCEPTION(common::InvalidOptionException("Mix of varying read lengths is not supported. Found: " +
                    boost::lexical_cast<std::string>(flowcellInfo_.readLengths_.second) + " and " + boost::lexical_cast<std::string>(block.getLSeq()) +
                    " for read 2"));
            }
        }

        // TODO: parse bam header to gather information about flowcells rather than rely on being able to meet
        // all flowcells in the first data block
        // scan until the end of the buffer or if we did not collect the information we need
        return !lastBlock || (paired_ && (!flowcellInfo_.readLengths_.first || !flowcellInfo_.readLengths_.second));
    }
private:
    std::string parseFlowcellId(const bam::BamBlockHeader &block)
    {
        const char *readNameEnd = block.read_name + block.getReadNameLength();
        const char *colon = std::find(block.read_name, readNameEnd, ':');
        if (readNameEnd == colon)
        {
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(std::string("Unable to parse flowcell id from read name. ") + block.read_name));
        }

        return std::string(block.read_name, colon);
    }
};


BamFlowcellInfo BamFlowcell::parseBamFlowcellInfo(
    const BamPath &laneFilePath,
    const bool allowVariableReadLength,
    const bool allowMixedFlowcells)
{
    BamFlowcellInfo ret;

    if (!laneFilePath.path_.empty())
    {
        common::ThreadVector threads(1);
        io::BamLoader bamLoader(0, threads, 1);
        bamLoader.open(laneFilePath.path_);

        MetadataParser metadataParser(ret);
        bamLoader.load
        (
            boost::make_tuple(
                boost::bind(&MetadataParser::parseMetadata, &metadataParser, allowVariableReadLength, allowMixedFlowcells, _1, _2),
                boost::bind(&nothing))
        );
    }

    ret.lanes_.push_back(laneFilePath.lane_);

    return ret;
}

flowcell::Layout BamFlowcell::createFilteredFlowcell(
    const bool detectSimpleIndels,
    const std::string &tilesFilter,
    const boost::filesystem::path &baseCallsDirectory,
    const unsigned laneNumberMax,
    std::string useBasesMask,
    const bool allowVariableReadLength,
    const std::string &seedDescriptor,
    const unsigned seedLength,
    const reference::ReferenceMetadataList &referenceMetadataList,
    unsigned &firstPassSeeds)
{

    BamPath flowcellFilePath = findBamPath(baseCallsDirectory);

    BamFlowcellInfo flowcellInfo = parseBamFlowcellInfo(flowcellFilePath, allowVariableReadLength,
                                                        "default" != useBasesMask && std::string::npos == useBasesMask.find('*'));

    std::vector<unsigned int> readLengths;
    if (flowcellInfo.readLengths_.first)
    {
        readLengths.push_back(flowcellInfo.readLengths_.first);
    }
    if (flowcellInfo.readLengths_.second)
    {
        readLengths.push_back(flowcellInfo.readLengths_.second);
    }

    if ("default" == useBasesMask)
    {
        if (readLengths.size() == 1)
        {
            useBasesMask = "y*";
        }
        else if (readLengths.size() == 2)
        {
            useBasesMask = "y*,y*";
        }
        else
        {
            const boost::format message =
                boost::format("\n   *** Could not guess the use-bases-mask for '%s', please supply the explicit value ***\n") %
                baseCallsDirectory.string();
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }
    }

    std::vector<unsigned int> readFirstCycles;

    ParsedUseBasesMask parsedUseBasesMask;
    alignment::SeedMetadataList seedMetadataList;
    if (!readLengths.empty())
    {
        parsedUseBasesMask = parseUseBasesMask(readFirstCycles, readLengths, seedLength, useBasesMask, baseCallsDirectory);
        seedMetadataList = parseSeedDescriptor(detectSimpleIndels, parsedUseBasesMask.dataReads_, seedDescriptor, seedLength, firstPassSeeds);
    }

    flowcell::Layout fc(baseCallsDirectory,
                        flowcell::Layout::Bam,
                        flowcell::BamFlowcellData(),
                        laneNumberMax,
                        std::vector<unsigned>(),
                        parsedUseBasesMask.dataReads_,
                        seedMetadataList, flowcellInfo.flowcellId_);

    std::string regexString(tilesFilter);
    std::replace(regexString.begin(), regexString.end(), ',', '|');
    boost::regex re(regexString);
    BOOST_FOREACH(const unsigned int lane, flowcellInfo.getLanes())
    {
        std::string laneString((boost::format("s_%d") % lane).str());
        if (boost::regex_search(laneString, re))
        {
            fc.addTile(lane, 1);
        }
    }

    return fc;
}


} // namespace alignOptions
} // namespace option
} // namespace isaac
