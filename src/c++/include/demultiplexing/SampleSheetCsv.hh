/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** BSD 2-Clause License
 **
 ** You should have received a copy of the BSD 2-Clause License
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
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
