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
 ** \file BclFlowcell.hh
 **
 ** Generate flowcell object out of BaseCalls/config.xml
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>

#include <boost/algorithm/string/regex.hpp>

#include "basecalls/ConfigXml.hh"

#include "alignOptions/UseBasesMaskOption.hh"
#include "io/FastqReader.hh"

#include "FastqFlowcell.hh"
#include "SeedDescriptorOption.hh"


namespace isaac
{
namespace options
{
namespace alignOptions
{

class CasavaFastqParser : boost::noncopyable
{
    const io::FastqReader& fastq_;
public:
    CasavaFastqParser(const io::FastqReader& fastq) :
        fastq_(fastq){}
    std::string parseFlowcellId()
    {
        ISAAC_ASSERT_MSG(fastq_.hasData(), "Can't parse empty fastq");
        io::FastqReader::IteratorPair header = fastq_.getHeader();
        static const std::string fastqHeaderDelimiters(": ");

        if ('@' != *header.first)
        {
            BOOST_THROW_EXCEPTION(io::FastqFormatException((boost::format("Fastq header must begin with @: %s") %
                std::string(header.first, header.second)).str()));
        }
        // separator between instrument name and run number
        io::FastqReader::IteratorPair::first_type flowcellIdBegin =
            std::find_first_of(header.first + 1, header.second,
                               fastqHeaderDelimiters.begin(), fastqHeaderDelimiters.end());

        if (header.second != flowcellIdBegin)
        {
            // separator between run number and flowcell id
            flowcellIdBegin = std::find_first_of(flowcellIdBegin + 1, header.second,
                                    fastqHeaderDelimiters.begin(), fastqHeaderDelimiters.end());

            if (header.second != flowcellIdBegin)
            {
                // separator between flowcell id and lane number
                const io::FastqReader::IteratorPair::first_type end =
                    std::find_first_of(flowcellIdBegin + 1, header.second,
                                       fastqHeaderDelimiters.begin(), fastqHeaderDelimiters.end());

                return std::string(flowcellIdBegin + 1, end);
            }
        }
        return "";
    }
    unsigned parseReadLength()
    {
        ISAAC_ASSERT_MSG(fastq_.hasData(), "Can't parse empty fastq");
        return fastq_.getReadLength();
    }
};


FastqFlowcell::FastqPathPairList FastqFlowcell::findFastqPathPairs(
    const flowcell::Layout::Format format,
    const boost::filesystem::path &baseCallsDirectory)
{
    FastqPathPairList ret;

    for (unsigned lane = 1; MAX_LANE_NUMBER >= lane; ++lane)
    {
        boost::filesystem::path r1Path;
        flowcell::Layout::getFastqFilePath(1, lane, baseCallsDirectory, flowcell::Layout::FastqGz == format, r1Path);

        FastqPathPair p;
        if (boost::filesystem::exists(r1Path))
        {
            p.r1Path_ = r1Path;
        }

        boost::filesystem::path r2Path;
        flowcell::Layout::getFastqFilePath(2, lane, baseCallsDirectory, flowcell::Layout::FastqGz == format, r2Path);

        if (boost::filesystem::exists(r2Path))
        {
            p.r2Path_ = r2Path;
        }

        if (!p.r1Path_.empty() || !p.r2Path_.empty())
        {
            p.lane_ = lane;
            ret.push_back(p);
        }
    }

    return ret;
}

FastqFlowcellInfo FastqFlowcell::parseFastqFlowcellInfo(
    const FastqPathPair &laneFilePaths)
{
    FastqFlowcellInfo ret;

    if (!laneFilePaths.r1Path_.empty())
    {
        io::FastqReader reader(false, laneFilePaths.r1Path_);
        CasavaFastqParser parser(reader);
        ret.readLengths_.first = parser.parseReadLength();

        ret.flowcellId_ = parser.parseFlowcellId();
    }
    if (!laneFilePaths.r2Path_.empty())
    {
        io::FastqReader reader(false, laneFilePaths.r2Path_);
        CasavaFastqParser parser(reader);
        ret.readLengths_.second = parser.parseReadLength();

        if (ret.flowcellId_.empty())
        {
            ret.flowcellId_ = parser.parseFlowcellId();
        }
        else if (ret.flowcellId_ != parser.parseFlowcellId())
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Flowcell ID mismatch between fastq reads %s vs %s, %s, %s") %
                ret.flowcellId_ % parser.parseFlowcellId() % laneFilePaths.r1Path_ % laneFilePaths.r2Path_).str()));
        }
    }

    ret.lanes_.push_back(laneFilePaths.lane_);

    return ret;
}

FastqFlowcellInfo FastqFlowcell::parseFastqFlowcellInfo(
    const FastqPathPairList &flowcellFilePaths)
{
    FastqFlowcellInfo ret = parseFastqFlowcellInfo(flowcellFilePaths.front());
    for (FastqPathPairList::const_iterator it = flowcellFilePaths.begin() + 1;
        flowcellFilePaths.end() != it; ++it)
    {
        FastqFlowcellInfo anotherLane = parseFastqFlowcellInfo(*it);
        if (anotherLane.flowcellId_ != ret.flowcellId_)
        {
            ISAAC_THREAD_CERR << "WARNING: Flowcell id mismatch across the lanes of the same flowcell" <<
                anotherLane << " vs " << ret << std::endl;
        }

        if (anotherLane.readLengths_ != ret.readLengths_)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Read lengths mismatch between lanes of the same flowcell %s vs %s") %
                anotherLane % ret).str()));
        }

        ret.lanes_.push_back(it->lane_);
    }
    ISAAC_THREAD_CERR << ret << std::endl;
    return ret;
}

flowcell::Layout FastqFlowcell::createFilteredFlowcell(
    const std::string &tilesFilter,
    const boost::filesystem::path &baseCallsDirectory,
    const flowcell::Layout::Format format,
    std::string useBasesMask,
    const std::string &seedDescriptor,
    const reference::ReferenceMetadataList &referenceMetadataList)
{

    FastqPathPairList flowcellFilePaths = findFastqPathPairs(format, baseCallsDirectory);
    if (flowcellFilePaths.empty())
    {
        const boost::format message = boost::format("\n   *** Could not find any fastq lanes in: %s ***\n") %
            baseCallsDirectory;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    FastqFlowcellInfo flowcellInfo = parseFastqFlowcellInfo(flowcellFilePaths);

    std::vector<unsigned int> readLengths;
    readLengths.push_back(flowcellInfo.readLengths_.first);
    readLengths.push_back(flowcellInfo.readLengths_.second);

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
        else
        {
            const boost::format message =
                boost::format("\n   *** Could not guess the use-bases-mask for '%s', please supply the explicit value ***\n") %
                baseCallsDirectory.string();
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }
    }

    std::vector<unsigned int> readFirstCycles;

    ParsedUseBasesMask parsedUseBasesMask = parseUseBasesMask(readFirstCycles, readLengths, useBasesMask, baseCallsDirectory);

    const alignment::SeedMetadataList seedMetadataList =
        parseSeedDescriptor(parsedUseBasesMask.dataReads_, seedDescriptor);

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

    if (fc.getLaneIds().empty())
    {
        const boost::format message = boost::format("\n   *** Could not find any lanes matching the '%s' in: %s ***\n") %
            tilesFilter % baseCallsDirectory;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    return fc;
}


} // namespace alignOptions
} // namespace option
} // namespace isaac
