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
 ** \file SampleSheetCsv.hh
 **
 ** SampleSheet.csv parser declaration
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_SAMPLE_SHEET_CSV_HH
#define iSAAC_DEMULTIPLEXING_SAMPLE_SHEET_CSV_HH

#include <boost/filesystem.hpp>

#include "basecalls/ConfigXml.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "reference/ReferenceMetadata.hh"

namespace isaac
{
namespace demultiplexing
{

flowcell::BarcodeMetadataList loadSampleSheetCsv(
    const boost::filesystem::path &sampleSheetPath,
    const std::string &assumedFlowcellId,
    const std::string &expectedFlowcellId,
    const unsigned expectedBarcodeLength,
    const unsigned flowcellIndex,
    const std::vector<unsigned> &lanes,
    const reference::ReferenceMetadataList &referenceMetadataList,
    const flowcell::SequencingAdapterMetadataList &defaultAdapters);

} // namespace demultiplexing
} // namespace isaac

#endif // #ifndef iSAAC_DEMULTIPLEXING_SAMPLE_SHEET_CSV_HH
