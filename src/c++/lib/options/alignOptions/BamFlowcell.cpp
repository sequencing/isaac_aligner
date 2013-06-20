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

BamFlowcell::BamPathList BamFlowcell::findBamPaths(
    const boost::filesystem::path &baseCallsDirectory)
{
    BamPathList ret;

    boost::filesystem::path path;
    flowcell::Layout::getBamFilePath(baseCallsDirectory, path);

    BamPath p;
    if (boost::filesystem::exists(path))
    {
        if (boost::filesystem::is_regular_file(path))
        {
            p.path_ = path;
            p.lane_ = 1;
            ret.push_back(p);
        }
        else
        {
            BOOST_THROW_EXCEPTION(common::InvalidOptionException("Bam --base-calls must be a regular file. Got: "  + path.string()));
        }
    }
    else
    {
        BOOST_THROW_EXCEPTION(common::InvalidOptionException("Bam file does not exist: "  + path.string()));
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

    bool parseMetadata(const bam::BamBlockHeader &block)
    {
        const std::string flowcellId = parseFlowcellId(block);
        if (flowcellInfo_.flowcellId_.empty())
        {
            flowcellInfo_.flowcellId_ = flowcellId;
        }
        else if (flowcellInfo_.flowcellId_ != flowcellId)
        {
            BOOST_THROW_EXCEPTION(common::InvalidOptionException("Multiple flowcell detected in the bam file: " +
                flowcellInfo_.flowcellId_ + " and " + flowcellId + ". Multiple flowcell data is not supported."));
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
            else if (flowcellInfo_.readLengths_.first != unsigned(block.getLSeq()))
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
            else if (flowcellInfo_.readLengths_.second != unsigned(block.getLSeq()))
            {
                BOOST_THROW_EXCEPTION(common::InvalidOptionException("Mix of varying read lengths is not supported. Found: " +
                    boost::lexical_cast<std::string>(flowcellInfo_.readLengths_.second) + " and " + boost::lexical_cast<std::string>(block.getLSeq()) +
                    " for read 2"));
            }
        }

        return paired_ && (!flowcellInfo_.readLengths_.first || !flowcellInfo_.readLengths_.second);
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
    const BamPath &laneFilePaths)
{
    BamFlowcellInfo ret;

    if (!laneFilePaths.path_.empty())
    {
        common::ThreadVector threads(1);
        io::BamLoader bamLoader(0, threads, 1);
        bamLoader.open(laneFilePaths.path_);

        MetadataParser metadataParser(ret);
        bamLoader.load
        (
            boost::make_tuple(
                boost::bind(&MetadataParser::parseMetadata, &metadataParser, _1),
                boost::bind(&nothing))
        );
    }

    ret.lanes_.push_back(laneFilePaths.lane_);

    return ret;
}

BamFlowcellInfo BamFlowcell::parseBamFlowcellInfo(
    const BamPathList &flowcellFilePaths,
    const bool allowVariableFastqLength)
{
    BamFlowcellInfo ret;
    bool flowcellInfoReady = false;
    for (BamPathList::const_iterator it = flowcellFilePaths.begin();
        flowcellFilePaths.end() != it; ++it)
    {
        BamFlowcellInfo anotherLane = parseBamFlowcellInfo(*it);
        if (!flowcellInfoReady)
        {
            if (anotherLane.readLengths_.first || anotherLane.readLengths_.second)
            {
                ret = anotherLane;
                flowcellInfoReady = true;
            }
            else
            {
                ISAAC_THREAD_CERR << "WARNING: Skipping lane " << it->lane_ << " due to read length 0" << std::endl;
            }
        }
        else
        {
            if (!anotherLane.readLengths_.first && !anotherLane.readLengths_.second)
            {
                ISAAC_THREAD_CERR << "WARNING: Skipping lane " << it->lane_ << " due to read length 0" << std::endl;
            }
            else
            {
                // With allowVariableFastqLength the read lengths are normally forced by use-bases-mask
                // Ignore discrepancy here.
                if (!allowVariableFastqLength && anotherLane.readLengths_ != ret.readLengths_)
                {
                    BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Read lengths mismatch between lanes of the same flowcell %s vs %s") %
                        anotherLane % ret).str()));
                }

                if (anotherLane.flowcellId_ != ret.flowcellId_)
                {
                    ISAAC_THREAD_CERR << "WARNING: Flowcell id mismatch across the lanes of the same flowcell" <<
                        anotherLane << " vs " << ret << std::endl;
                }

                ret.lanes_.push_back(it->lane_);
            }
        }
    }
    ISAAC_THREAD_CERR << ret << std::endl;
    return ret;
}

flowcell::Layout BamFlowcell::createFilteredFlowcell(
    const bool detectSimpleIndels,
    const std::string &tilesFilter,
    const boost::filesystem::path &baseCallsDirectory,
    const flowcell::Layout::Format format,
    std::string useBasesMask,
    const bool allowVariableFastqLength,
    const std::string &seedDescriptor,
    const unsigned seedLength,
    const reference::ReferenceMetadataList &referenceMetadataList,
    unsigned &firstPassSeeds)
{

    ISAAC_ASSERT_MSG(flowcell::Layout::Bam == format, "Wrong format passed to BamFlowcell::createFilteredFlowcell");
    BamPathList flowcellFilePaths = findBamPaths(baseCallsDirectory);
    ISAAC_ASSERT_MSG(!flowcellFilePaths.empty(), "Missing bam file should have caused an exception in findBamPaths");

    BamFlowcellInfo flowcellInfo = parseBamFlowcellInfo(flowcellFilePaths, allowVariableFastqLength);

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
                        format,
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
