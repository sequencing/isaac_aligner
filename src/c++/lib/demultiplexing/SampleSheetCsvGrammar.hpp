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
 ** \file SampleSheetCsvGrammar.cpp
 **
 ** SampleSheet.csv grammar definition
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_SAMPLE_SHEET_CSV_GRAMMAR_HH
#define iSAAC_DEMULTIPLEXING_SAMPLE_SHEET_CSV_GRAMMAR_HH

#include "StandardSampleSheetCsvGrammar.hpp"
#include "MiSeqDualBarcodeSampleSheetCsvGrammar.hpp"
#include "MiSeqSingleBarcodeSampleSheetCsvGrammar.hpp"
#include "MiSeqNonMultiplexedSampleSheetCsvGrammar.hpp"

namespace isaac
{
namespace demultiplexing
{

namespace bs=boost::spirit;

template <typename Iterator>
struct SampleSheetCsvGrammar :
    bs::qi::grammar<Iterator, flowcell::BarcodeMetadataList()>
{
    SampleSheetCsvGrammar(const flowcell::SequencingAdapterMetadataList &defaultAdapters) :
        SampleSheetCsvGrammar::base_type(start_),
        miSeqDualBarcodeGrammar_(defaultAdapters),
        miSeqSingleBarcodeGrammar_(defaultAdapters),
        miSeqNonMultiplexedGrammar_(defaultAdapters),
        standardGrammar_(defaultAdapters)
    {
        start_ = miSeqDualBarcodeGrammar_.start_ | miSeqSingleBarcodeGrammar_.start_ |
            miSeqNonMultiplexedGrammar_.start_ | standardGrammar_.start_;
    }

    MiSeqDualBarcodeSampleSheetCsvGrammar<Iterator> miSeqDualBarcodeGrammar_;
    MiSeqSingleBarcodeSampleSheetCsvGrammar<Iterator> miSeqSingleBarcodeGrammar_;
    MiSeqNonMultiplexedSampleSheetCsvGrammar<Iterator> miSeqNonMultiplexedGrammar_;
    StandardSampleSheetCsvGrammar<Iterator> standardGrammar_;

    bs::qi::rule<Iterator, flowcell::BarcodeMetadataList()> start_;
};

} // namespace demultiplexing
} // namespace isaac

#endif //iSAAC_DEMULTIPLEXING_SAMPLE_SHEET_CSV_GRAMMAR_HH
