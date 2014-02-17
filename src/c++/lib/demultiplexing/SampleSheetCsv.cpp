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
 ** \file SampleSheetCsv.cpp
 **
 ** SampleSheet.csv parser implementation
 **
 ** \author Roman Petrovski
 **/
#include <boost/foreach.hpp>
#include "demultiplexing/SampleSheetCsv.hh"
#include "SampleSheetCsvGrammar.hpp"

namespace isaac
{
namespace demultiplexing
{

void loadSampleSheetCsv(
    const boost::filesystem::path &sampleSheetPath,
    const flowcell::SequencingAdapterMetadataList &defaultAdapters,
    flowcell::BarcodeMetadataList &sampleSheet)
{
    std::ifstream is(sampleSheetPath.c_str());
    if (!is) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open sample sheet file " + sampleSheetPath.string()));
    }
    is.unsetf(std::ios::skipws);
    std::vector<char> sampleSheetBuffer;
    std::copy(std::istream_iterator<char>(is), std::istream_iterator<char>(), std::back_inserter(sampleSheetBuffer));
    if (!is.eof())
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to read " + sampleSheetPath.string() + " to the end"));
    }
    std::vector<char>::const_iterator parseIt(sampleSheetBuffer.begin());
    std::vector<char>::const_iterator parseEnd(sampleSheetBuffer.end());
    SampleSheetCsvGrammar<std::vector<char>::const_iterator> parser(defaultAdapters);

    try
    {
        if (!boost::spirit::qi::parse(parseIt, parseEnd, parser, sampleSheet) || parseEnd != parseIt)
        {
            const std::string message = "Could not parse the sample sheet csv stream text:\n" +
                std::string(parseIt, parseEnd);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message));
        }
    }
    catch(boost::spirit::qi::expectation_failure<std::vector<char>::const_iterator> const &e)
    {
        const std::vector<char>::const_iterator bufferBegin(sampleSheetBuffer.begin());

        const std::string message = (boost::format("Could not parse the sample sheet csv. "
            "Expected:"
            "\n==================================================================\n"
            "%s"
            "\n==================================================================\n"
            "Got:"
            "\n==================================================================\n"
            "%s"
            "\n==================================================================\n"
            " at offset: %d") % e.what_ % std::string(e.first, e.last) %
            std::distance(bufferBegin, e.first)).str();
        BOOST_THROW_EXCEPTION(common::IoException(errno, message));
    }

    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, sampleSheet)
    {
        ISAAC_THREAD_CERR << barcode << std::endl;
    }
}

/**
 * \param sampleSheetPath If empty, the default non-indexed barcode will be generated for each lane present in the
 *                        config.xml
 * \param expectedBarcodeLength barcode length inferred from base calls metadata and command line
 */
flowcell::BarcodeMetadataList loadSampleSheetCsv(
    const boost::filesystem::path &sampleSheetPath,
    const std::string &assumedFlowcellId,
    const std::string &expectedFlowcellId,
    const unsigned expectedBarcodeLength,
    const unsigned flowcellIndex,
    const std::vector<unsigned> &lanes,
    const reference::ReferenceMetadataList &referenceMetadataList,
    const flowcell::SequencingAdapterMetadataList &defaultAdapters)
{
    flowcell::BarcodeMetadataList ret;
    const reference::ReferenceMetadataList::const_iterator defaultReference =
        std::find_if(referenceMetadataList.begin(), referenceMetadataList.end(),
                     boost::bind(&reference::ReferenceMetadata::getName, _1) == "default");

    if (sampleSheetPath.empty())
    {
        ISAAC_ASSERT_MSG(referenceMetadataList.end() != defaultReference,
                         "When no sample sheet is provided, 'default' reference must have a specification");
        // push the 'no index' barcode for each lane
        BOOST_FOREACH(const unsigned lane, lanes)
        {
            ret.push_back(flowcell::BarcodeMetadata::constructNoIndexBarcode(
                assumedFlowcellId, flowcellIndex, lane, defaultReference->getIndex(), defaultAdapters));
            ISAAC_THREAD_CERR << "Generated 'none' barcode: " << ret.back() << std::endl;

        }
        return ret;
    }

    loadSampleSheetCsv(sampleSheetPath, defaultAdapters, ret);
    if (ret.empty())
    {
        const boost::format message =
            boost::format("\n   *** Sample sheet must define at least one sample. "
                "Sample sheet: %s***\n") % sampleSheetPath.string();
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    std::vector<unsigned> lanesWithoutUnknownBarcodeSpec(lanes);
    BOOST_FOREACH(flowcell::BarcodeMetadata &barcode, ret)
    {
        // TODO: WARNING if the '-' in the barcode sequence is not at the read border
        // WARN if flowcell id is available from config.xml and does not match
        if (!expectedFlowcellId.empty() && barcode.getFlowcellId() != expectedFlowcellId)
        {
            ISAAC_THREAD_CERR << "WARNING: the flowcell id " << expectedFlowcellId << " from BaseCalls data"
                << " does not match the one read for the barcode " << barcode << " from " << sampleSheetPath << std::endl;

            // Make sure the flowcell id matches or else there will be broken xmls somewhere
        }
        if (!barcode.isUnknown() && barcode.getSequenceLength() != expectedBarcodeLength)
        {
            const boost::format message =
                boost::format("\n   *** Sample sheet barcode sequence length does not match barcode cycles data. "
                    "Expected %d, got: %s. Sample sheet: %s***\n") %
                expectedBarcodeLength % barcode % sampleSheetPath.string();
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }

        // force assumed flowcell id to make sure there is consistency between flowcell metadata and barcode metadata
        // in the system
        barcode.setFlowcellId(assumedFlowcellId);
        barcode.setFlowcellIndex(flowcellIndex);
        if (barcode.isDefault())
        {
            lanesWithoutUnknownBarcodeSpec.erase(
                std::remove(lanesWithoutUnknownBarcodeSpec.begin(), lanesWithoutUnknownBarcodeSpec.end(), barcode.getLane()),
                lanesWithoutUnknownBarcodeSpec.end());
        }

        reference::ReferenceMetadataList::const_iterator barcodeReference =
            std::find_if(referenceMetadataList.begin(), referenceMetadataList.end(),
                         boost::bind(&reference::ReferenceMetadata::getName, _1) == barcode.getReference());

        if (referenceMetadataList.end() != barcodeReference)
        {
            barcode.setReferenceIndex(barcodeReference->getIndex());
        }
        else if (referenceMetadataList.end() != defaultReference)
        {
            barcode.setReferenceIndex(defaultReference->getIndex());
        }
    }

    const reference::ReferenceMetadataList::const_iterator unknownReference =
        std::find_if(referenceMetadataList.begin(), referenceMetadataList.end(),
                     boost::bind(&reference::ReferenceMetadata::getName, _1) == "unknown");

    // push the default 'unknown indexes' barcode for each lane that does not have it specified in sample sheet.
    BOOST_FOREACH(const unsigned lane, lanesWithoutUnknownBarcodeSpec)
    {
        ret.push_back(flowcell::BarcodeMetadata::constructUnknownBarcode(
            assumedFlowcellId, flowcellIndex, lane,
            referenceMetadataList.end() == unknownReference
                ? flowcell::BarcodeMetadata::UNMAPPED_REFERENCE_INDEX : unknownReference->getIndex(), defaultAdapters));

        ISAAC_THREAD_CERR << "Generated 'unknown' barcode: " << ret.back() << std::endl;
    }

    return ret;
}

} // namespace demultiplexing
} // namespace isaac

