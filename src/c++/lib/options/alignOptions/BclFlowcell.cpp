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
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/support/unused.hpp>

#include <boost/algorithm/string/regex.hpp>

#include "basecalls/ConfigXml.hh"

#include "alignOptions/UseBasesMaskOption.hh"

#include "BclFlowcell.hh"
#include "SeedDescriptorOption.hh"


namespace isaac
{
namespace options
{
namespace alignOptions
{

static basecalls::ConfigXml parseBasecallsConfigXml(const boost::filesystem::path &baseCallsDirectory)
{
    const boost::filesystem::path basecallsConfigXml(baseCallsDirectory/"config.xml");
    std::ifstream is(basecallsConfigXml.string().c_str());
    if (!is) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open basecalls config file " + basecallsConfigXml.string()));
    }
    basecalls::ConfigXml cfg;
    if (!(is >> cfg)) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to read from basecalls config file " + basecallsConfigXml.string()));
    }
    return cfg;
}

flowcell::Layout BclFlowcell::createFilteredFlowcell(
    const std::string &tilesFilter,
    const boost::filesystem::path &baseCallsDirectory,
    const flowcell::Layout::Format format,
    std::string useBasesMask,
    const std::string &seedDescriptor,
    const reference::ReferenceMetadataList &referenceMetadataList)
{
    basecalls::ConfigXml cfg = parseBasecallsConfigXml(baseCallsDirectory);

    using boost::phoenix::at;
    using boost::phoenix::bind;
    using boost::phoenix::arg_names::_1;
    using boost::phoenix::ref;
    std::vector<basecalls::ConfigXml::RunParametersRead> cfgReads(cfg.getRunParametersReads());

    std::vector<unsigned int> readLengths;
    std::transform(cfgReads.begin(), cfgReads.end(),
                   std::back_inserter(readLengths),
                   bind(&basecalls::ConfigXml::RunParametersRead::lastCycle_, _1) -
                       bind(&basecalls::ConfigXml::RunParametersRead::firstCycle_, _1) + 1);

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
                baseCallsDirectory.string();
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }

    }

    std::vector<unsigned int> readFirstCycles;
    std::transform(cfgReads.begin(), cfgReads.end(),
                   std::back_inserter(readFirstCycles),
                   bind(&basecalls::ConfigXml::RunParametersRead::firstCycle_, _1));


    ParsedUseBasesMask parsedUseBasesMask = parseUseBasesMask(readFirstCycles, readLengths, useBasesMask, baseCallsDirectory);
    std::vector<unsigned> barcodeCycles;
    BOOST_FOREACH(const flowcell::ReadMetadata &barcodeRead, parsedUseBasesMask.indexReads_)
    {
        barcodeCycles.insert(barcodeCycles.end(), barcodeRead.getCycles().begin(), barcodeRead.getCycles().end());
    }

    const alignment::SeedMetadataList seedMetadataList =
        parseSeedDescriptor(parsedUseBasesMask.dataReads_, seedDescriptor);

    flowcell::Layout fc(baseCallsDirectory,
                        format,
                        barcodeCycles,
                        parsedUseBasesMask.dataReads_,
                        seedMetadataList, cfg.getFlowcellId());

    std::string regexString(tilesFilter);
    std::replace(regexString.begin(), regexString.end(), ',', '|');
    boost::regex re(regexString);
    BOOST_FOREACH(const unsigned int lane, cfg.getLanes())
    {
        BOOST_FOREACH(const unsigned int tile, cfg.getTiles(lane))
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
            tilesFilter % baseCallsDirectory;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    return fc;
}


} // namespace alignOptions
} // namespace option
} // namespace isaac
