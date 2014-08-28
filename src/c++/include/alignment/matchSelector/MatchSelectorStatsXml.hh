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
 ** \file MatchSelectorStatsXml.hh
 **
 ** \brief Xml Serialization of MatchSelector statistics.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_ALIGNMENT_MATCH_SELECTOR_STATS_XML_H
#define ISAAC_ALIGNMENT_MATCH_SELECTOR_STATS_XML_H

#include "MatchSelectorStats.hh"
#include "xml/XmlWriter.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

class MatchSelectorStatsXml
{
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::TileMetadataList &tileMetadataList_;
    const std::vector<MatchSelectorStats> &stats_;
public:
    MatchSelectorStatsXml(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::TileMetadataList &tileMetadataList,
        const std::vector<MatchSelectorStats> &stats) :
            flowcellLayoutList_(flowcellLayoutList),
            barcodeMetadataList_(barcodeMetadataList),
            tileMetadataList_(tileMetadataList),
            stats_(stats)
    {}

    void serialize(std::ostream& os) const;

private:
    void serializeBarcodes(
        xml::XmlWriter &xmlWriter,
        const flowcell::Layout &flowcell) const;

    void serializeReads(
        xml::XmlWriter &xmlWriter,
        const flowcell::Layout &flowcell) const;

    void serlializeTile(
        xml::XmlWriter &xmlWriter,
        const flowcell::TileMetadata &tile) const;

    void serlializeTileRead(
        xml::XmlWriter &xmlWriter,
        const flowcell::ReadMetadata &read,
        const TileStats& tileStats) const;

    void serlializeTemplateAlignmentScoreDistribution(
        xml::XmlWriter &xmlWriter,
        const TileStats& tileStats) const;

    void serializeTileBarcode(
        xml::XmlWriter &xmlWriter,
        const unsigned read,
        const bool lowestRead,
        const bool passesFilter,
        const TileBarcodeStats& tileStats) const;

    typedef std::map<unsigned, boost::array<matchSelector::TileBarcodeStats, 2> > ReadBarcodeStats;
    typedef std::map<unsigned, ReadBarcodeStats> TileReadBarcodeStats;
    typedef std::map<std::string, TileReadBarcodeStats> BarcodeTileReadBarcodeStats;
    typedef std::map<std::string, BarcodeTileReadBarcodeStats> SampleBarcodeTileReadBarcodeStats;
    typedef std::map<std::string, SampleBarcodeTileReadBarcodeStats> ProjectSampleBarcodeTileReadBarcodeStats;
    typedef std::map<std::string, ProjectSampleBarcodeTileReadBarcodeStats> FlowcellProjectSampleBarcodeTileReadBarcodeStats;

    void serializeFlowcellProjects(
        xml::XmlWriter &xmlWriter,
        const std::vector<unsigned> &allLanes,
        const ProjectSampleBarcodeTileReadBarcodeStats &flowcellProjectSampleBarcodeStats) const;
};

} //namespace matchSelector
} //namespace alignment
} //namespace isaac

#endif //ISAAC_ALIGNMENT_MATCH_SELECTOR_STATS_XML_H
