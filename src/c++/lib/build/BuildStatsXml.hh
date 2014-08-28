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
 ** \file BuildStatsXml.hh
 **
 ** \brief Xml Serialization of Build statistics.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_BUILD_BUILD_STATS_XML_H
#define ISAAC_BUILD_BUILD_STATS_XML_H

#include "alignment/BinMetadata.hh"
#include "build/BuildStats.hh"
#include "reference/SortedReferenceMetadata.hh"
#include "xml/XmlWriter.hh"

namespace isaac
{
namespace build
{

class BuildStatsXml
{
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList_;
    const alignment::BinMetadataCRefList &bins_;
    flowcell::BarcodeMetadataList orderedBarcodeMetadataList_;
    const BuildStats &buildStats_;

    void dumpContigs(
        xml::XmlWriter &xmlWriter,
        const std::vector<reference::SortedReferenceMetadata::Contig> &contigs,
        flowcell::BarcodeMetadataList::const_iterator sampleBarcodesBegin,
        flowcell::BarcodeMetadataList::const_iterator sampleBarcodesEnd);

public:
    BuildStatsXml(
        const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const alignment::BinMetadataCRefList &bins,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const BuildStats &buildStats);

    void serialize(std::ostream &os);
};

} //namespace build
} //namespace isaac

#endif //ISAAC_BUILD_BUILD_STATS_XML_H
