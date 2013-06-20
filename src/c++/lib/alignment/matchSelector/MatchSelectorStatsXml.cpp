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
 ** \file MatchSelectorStatsXml.cpp
 **
 ** \brief Xml Serialization of MatchSelector statistics.
 **
 ** \author Roman Petrovski
 **/

#include <boost/lexical_cast.hpp>

#include "alignment/matchSelector/MatchSelectorStatsXml.hh"


namespace isaac
{
namespace alignment
{
namespace matchSelector
{

const std::string &alignmentModelName(const TemplateLengthStatistics::AlignmentModel alignmentModel)
{
    static const std::vector<std::string> modelNames = boost::assign::list_of
        ("FFp")("FRp")("RFp")("RRp")("FFm")("FRm")("RFm")("RRm")("unknown");
    ISAAC_ASSERT_MSG(modelNames.size() > unsigned(alignmentModel), "Incorrect alignmentModel value");

    return modelNames[alignmentModel];
}

const std::string &alignmentClassName(const TemplateLengthStatistics::AlignmentModel alignmentModel)
{
    static const std::vector<std::string> classNames = boost::assign::list_of
        ("Fp")("Rp")("Rm")("Fm")("Fm")("Rm")("Rp")("Fp")("unknown");

    ISAAC_ASSERT_MSG(classNames.size() > unsigned(alignmentModel), "Incorrect alignmentClass value");

    return classNames[alignmentModel];
}

void MatchSelectorStatsXml::serializeBarcodes(
    xml::XmlWriter &xmlWriter, const flowcell::Layout &flowcell) const
{
    std::string lastBarcodeName;
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList_)
    {
        if (barcode.getFlowcellIndex() == flowcell.getIndex())
        {
            if (barcode.getName() != lastBarcodeName)
            {
                if (!lastBarcodeName.empty())
                {
                    xmlWriter.endElement();
                }
                lastBarcodeName = barcode.getName();
                xmlWriter.startElement("Barcode");
                xmlWriter.writeAttribute("name", lastBarcodeName);
            }

            ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Lane")
            {
                xmlWriter.writeAttribute("number", barcode.getLane());
                xmlWriter.writeElement("ReferenceName", barcode.getReference());
            }
        }
    }
    if (!lastBarcodeName.empty())
    {
        xmlWriter.endElement();
    }
}

void MatchSelectorStatsXml::serializeReads(
    xml::XmlWriter &xmlWriter, const flowcell::Layout &flowcell) const
{

    BOOST_FOREACH(const flowcell::ReadMetadata& read, flowcell.getReadMetadataList())
    {
        ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Read")
        {
            xmlWriter.writeAttribute("number", read.getNumber());
            xmlWriter.writeElement("Length", read.getLength());
        }
    }

}

void MatchSelectorStatsXml::serializeFlowcellProjects(
    xml::XmlWriter &xmlWriter,
    const std::vector<unsigned> &allLanes,
    const ProjectSampleBarcodeTileReadBarcodeStats &flowcellProjectSampleBarcodeStats) const
{
    BOOST_FOREACH(const ProjectSampleBarcodeTileReadBarcodeStats::value_type &projectStats, flowcellProjectSampleBarcodeStats)
    {
        ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Project")
        {
            xmlWriter.writeAttribute("name", projectStats.first);
            BOOST_FOREACH(const SampleBarcodeTileReadBarcodeStats::value_type &sampleStats, projectStats.second)
            {
                ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Sample")
                {
                    xmlWriter.writeAttribute("name", sampleStats.first);
                    BOOST_FOREACH(const BarcodeTileReadBarcodeStats::value_type &barcodeStats, sampleStats.second)
                    {
                        ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Barcode")
                        {
                            xmlWriter.writeAttribute("name", barcodeStats.first);
                            unsigned lastCreatedLane = 0;
                            BOOST_FOREACH(const unsigned lane, allLanes)
                            {
                                BOOST_FOREACH(const TileReadBarcodeStats::value_type &tileStats, barcodeStats.second)
                                {
                                    if (tileMetadataList_.at(tileStats.first).getLane() == lane)
                                    {
                                        bool dumpReadIndependentStats = true;
                                        BOOST_FOREACH(const ReadBarcodeStats::value_type &readStats, tileStats.second)
                                        {
                                            if (lastCreatedLane != lane)
                                            {
                                                // avoid creating empty Lane elements as they screw up html reports
                                                xmlWriter.startElement("Lane");
                                                xmlWriter.writeAttribute("number", lane);
                                                lastCreatedLane = lane;
                                            }
                                            ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Tile")
                                            {
                                                xmlWriter.writeAttribute("number", tileMetadataList_.at(tileStats.first).getTile());
                                                serializeTileBarcode(xmlWriter, readStats.first, dumpReadIndependentStats, true, readStats.second[0]);
                                                serializeTileBarcode(xmlWriter, readStats.first, dumpReadIndependentStats, false, readStats.second[1]);
                                                dumpReadIndependentStats = false;
                                            }
                                        }
                                    }
                                }
                                if (lastCreatedLane == lane)
                                {
                                    xmlWriter.endElement();
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void MatchSelectorStatsXml::serialize(std::ostream& os) const
{
    xml::XmlWriter xmlWriter(os);

    FlowcellProjectSampleBarcodeTileReadBarcodeStats flowcellProjectSampleBarcodeStats;
    std::vector<unsigned> allLanes;

    ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Stats")
    {
        BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList_)
        {
            std::vector<unsigned int> lanes = flowcell.getLaneIds();
            allLanes.insert(allLanes.end(), lanes.begin(), lanes.end());
            ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Flowcell")
            {
                xmlWriter.writeAttribute("flowcell-id", flowcell.getFlowcellId());

                serializeReads(xmlWriter, flowcell);
                serializeBarcodes(xmlWriter, flowcell);

                unsigned lastCreatedLane = 0;
                BOOST_FOREACH(const unsigned lane, flowcell.getLaneIds())
                {
                    BOOST_FOREACH(const flowcell::TileMetadata& tile, tileMetadataList_)
                    {
                        if (lane == tile.getLane() && flowcell.getIndex() == tile.getFlowcellIndex())
                        {
                            if (lastCreatedLane != lane)
                            {
                                // avoid creating empty Lane elements as they screw up html reports
                                xmlWriter.startElement("Lane");
                                xmlWriter.writeAttribute("number", lane);
                                lastCreatedLane = lane;
                            }
                            serlializeTile(xmlWriter, tile);
                        }
                    }
                    if (lastCreatedLane == lane)
                    {
                        xmlWriter.endElement();
                    }
                }

                BOOST_FOREACH(const flowcell::TileMetadata& tile, tileMetadataList_)
                {
                    if (flowcell.getIndex() == tile.getFlowcellIndex())
                    {
                        BOOST_FOREACH(const flowcell::ReadMetadata& read, flowcellLayoutList_.at(tile.getFlowcellIndex()).getReadMetadataList())
                        {
                            BOOST_FOREACH(const flowcell::BarcodeMetadata& barcode, barcodeMetadataList_)
                            {
                                if (barcode.getFlowcellId() == tile.getFlowcellId() && barcode.getLane() == tile.getLane())
                                {
                                    const matchSelector::TileBarcodeStats &pfStat = stats_.at(tile.getIndex()).getReadBarcodeTileStat(read, barcode, true);
                                    flowcellProjectSampleBarcodeStats[barcode.getFlowcellId()][barcode.getProject()][barcode.getSampleName()][barcode.getName()][tile.getIndex()][read.getNumber()][0] += pfStat;
                                    flowcellProjectSampleBarcodeStats[barcode.getFlowcellId()][barcode.getProject()][barcode.getSampleName()]["all"][tile.getIndex()][read.getNumber()][0] += pfStat;
                                    flowcellProjectSampleBarcodeStats[barcode.getFlowcellId()][barcode.getProject()]["all"]["all"][tile.getIndex()][read.getNumber()][0] += pfStat;
                                    flowcellProjectSampleBarcodeStats[barcode.getFlowcellId()]["all"]["all"]["all"][tile.getIndex()][read.getNumber()][0] += pfStat;
                                    flowcellProjectSampleBarcodeStats["all"][barcode.getProject()][barcode.getSampleName()]["all"][tile.getIndex()][read.getNumber()][0] += pfStat;
                                    flowcellProjectSampleBarcodeStats["all"][barcode.getProject()]["all"]["all"][tile.getIndex()][read.getNumber()][0] += pfStat;
                                    flowcellProjectSampleBarcodeStats["all"]["all"]["all"]["all"][tile.getIndex()][read.getNumber()][0] += pfStat;

                                    const matchSelector::TileBarcodeStats &rawStat = stats_.at(tile.getIndex()).getReadBarcodeTileStat(read, barcode, false);
                                    flowcellProjectSampleBarcodeStats[barcode.getFlowcellId()][barcode.getProject()][barcode.getSampleName()][barcode.getName()][tile.getIndex()][read.getNumber()][1] += rawStat;
                                    flowcellProjectSampleBarcodeStats[barcode.getFlowcellId()][barcode.getProject()][barcode.getSampleName()]["all"][tile.getIndex()][read.getNumber()][1] += rawStat;
                                    flowcellProjectSampleBarcodeStats[barcode.getFlowcellId()][barcode.getProject()]["all"]["all"][tile.getIndex()][read.getNumber()][1] += rawStat;
                                    flowcellProjectSampleBarcodeStats[barcode.getFlowcellId()]["all"]["all"]["all"][tile.getIndex()][read.getNumber()][1] += rawStat;
                                    flowcellProjectSampleBarcodeStats["all"][barcode.getProject()][barcode.getSampleName()]["all"][tile.getIndex()][read.getNumber()][1] += rawStat;
                                    flowcellProjectSampleBarcodeStats["all"][barcode.getProject()]["all"]["all"][tile.getIndex()][read.getNumber()][1] += rawStat;
                                    flowcellProjectSampleBarcodeStats["all"]["all"]["all"]["all"][tile.getIndex()][read.getNumber()][1] += rawStat;
                                }
                            }
                        }
                    }
                }

                serializeFlowcellProjects(xmlWriter, flowcell.getLaneIds(), flowcellProjectSampleBarcodeStats.at(flowcell.getFlowcellId()));
            }
        }

        std::sort(allLanes.begin(), allLanes.end());
        allLanes.erase(std::unique(allLanes.begin(), allLanes.end()), allLanes.end());
        ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Flowcell")
        {
            xmlWriter.writeAttribute("flowcell-id", "all");
            serializeFlowcellProjects(xmlWriter, allLanes, flowcellProjectSampleBarcodeStats.at("all"));
        }
    }
}


void MatchSelectorStatsXml::serlializeTile(
    xml::XmlWriter &xmlWriter,
    const flowcell::TileMetadata &tile) const
{
    ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Tile")
    {
        xmlWriter.writeAttribute("number", tile.getTile());
        ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Pf")
        {
            BOOST_FOREACH(const flowcell::ReadMetadata& read, flowcellLayoutList_.at(tile.getFlowcellIndex()).getReadMetadataList())
            {
                serlializeTileRead(xmlWriter, read, stats_.at(tile.getIndex()).getReadTileStat(read, true));
            }
        }
        ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Raw")
        {
            BOOST_FOREACH(const flowcell::ReadMetadata& read, flowcellLayoutList_.at(tile.getFlowcellIndex()).getReadMetadataList())
            {
                serlializeTileRead(xmlWriter, read, stats_.at(tile.getIndex()).getReadTileStat(read, false));
            }
        }
    }
}

void MatchSelectorStatsXml::serlializeTileRead(
    xml::XmlWriter &xmlWriter,
    const flowcell::ReadMetadata &read,
    const TileStats& tileStats) const
{
    if (!read.getIndex())
    {
        serlializeTemplateAlignmentScoreDistribution(xmlWriter, tileStats);
    }

    ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Read")
    {
        xmlWriter.writeAttribute("number", read.getNumber());

        BOOST_STATIC_ASSERT(sizeof(tileStats.alignmentScoreMismatches_) == sizeof(tileStats.alignmentScoreFragments_));
        ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "AllFragments")
        {
            ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "AlignmentScoreDistribution")
            {
                for(unsigned long score = 0; score <= tileStats.maxAlignmentScore_; ++score)
                {
                    if (tileStats.alignmentScoreFragments_[score] || tileStats.alignmentScoreMismatches_[score])
                    {
                        ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Score")
                        {
                            xmlWriter.writeAttribute("fragments", tileStats.alignmentScoreFragments_[score]);
                            xmlWriter.writeAttribute("mismatches", tileStats.alignmentScoreMismatches_[score]);
                            xmlWriter.writeAttribute("score", score);
                        }
                    }
                }
            }
            ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "MismatchesByCycle")
            {
                for(unsigned long cycle = read.getFirstCycle(); cycle <= read.getLastCycle(); ++cycle)
                {
                    ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Cycle")
                    {
                        xmlWriter.writeAttribute("blanks", tileStats.cycleBlanks_[cycle]);
                        xmlWriter.writeAttribute("mismatches", tileStats.cycleMismatches_[cycle]);
                        xmlWriter.writeAttribute("number", cycle);
                    }
                }
            }

            ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "FragmentMismatchesByCycle")
            {
                for(unsigned long cycle = read.getFirstCycle(); cycle <= read.getLastCycle(); ++cycle)
                {
                    ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Cycle")
                    {
                        xmlWriter.writeAttribute("one", tileStats.cycle1MismatchFragments_[cycle]);
                        xmlWriter.writeAttribute("two", tileStats.cycle2MismatchFragments_[cycle]);
                        xmlWriter.writeAttribute("three", tileStats.cycle3MismatchFragments_[cycle]);
                        xmlWriter.writeAttribute("four", tileStats.cycle4MismatchFragments_[cycle]);
                        xmlWriter.writeAttribute("more", tileStats.cycleMoreMismatchFragments_[cycle]);
                        xmlWriter.writeAttribute("number", cycle);
                    }
                }
            }
        }

        ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "UniquelyAlignedFragments")
        {
            ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "MismatchesByCycle")
            {
                for(unsigned long cycle = read.getFirstCycle(); cycle <= read.getLastCycle(); ++cycle)
                {
                    ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Cycle")
                    {
                        xmlWriter.writeAttribute("blanks", tileStats.cycleUniquelyAlignedBlanks_[cycle]);
                        xmlWriter.writeAttribute("mismatches", tileStats.cycleUniquelyAlignedMismatches_[cycle]);
                        xmlWriter.writeAttribute("number", cycle);
                    }
                }
            }

            ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "FragmentMismatchesByCycle")
            {
                for(unsigned long cycle = read.getFirstCycle(); cycle <= read.getLastCycle(); ++cycle)
                {
                    ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Cycle")
                    {
                        xmlWriter.writeAttribute("one", tileStats.cycleUniquelyAligned1MismatchFragments_[cycle]);
                        xmlWriter.writeAttribute("two", tileStats.cycleUniquelyAligned2MismatchFragments_[cycle]);
                        xmlWriter.writeAttribute("three", tileStats.cycleUniquelyAligned3MismatchFragments_[cycle]);
                        xmlWriter.writeAttribute("four", tileStats.cycleUniquelyAligned4MismatchFragments_[cycle]);
                        xmlWriter.writeAttribute("more", tileStats.cycleUniquelyAlignedMoreMismatchFragments_[cycle]);
                        xmlWriter.writeAttribute("number", cycle);
                    }
                }
            }
            xmlWriter.writeElement("Count", tileStats.uniquelyAlignedFragmentCount_);
        }
    }
}

void MatchSelectorStatsXml::serlializeTemplateAlignmentScoreDistribution(
    xml::XmlWriter &xmlWriter,
    const TileStats& tileStats) const
{
    BOOST_STATIC_ASSERT(sizeof(tileStats.alignmentScoreTemplateMismatches_) == sizeof(tileStats.alignmentScoreTemplates_));
    ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "AlignmentScoreDistribution")
    {
        for(unsigned long score = 0; score <= tileStats.maxAlignmentScore_; ++score)
        {
            if (tileStats.alignmentScoreTemplates_[score] || tileStats.alignmentScoreMismatches_[score])
            {
                ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Score")
                {
                    xmlWriter.writeAttribute("templates", tileStats.alignmentScoreTemplates_[score]);
                    xmlWriter.writeAttribute("mismatches", tileStats.alignmentScoreTemplateMismatches_[score]);
                    xmlWriter.writeAttribute("score", score);
                }
            }
        }
    }
}

void MatchSelectorStatsXml::serializeTileBarcode(
    xml::XmlWriter &xmlWriter,
    const unsigned read,
    const bool dumpReadIndependentStats,
    const bool passesFilter,
    const TileBarcodeStats& tileStats) const
{
    // template length stats stored as read index 0 tiles with passesFilter == false
    if (dumpReadIndependentStats && !passesFilter)
    {
        ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "AssumedTemplateLength")
        {
            xmlWriter.writeElement("Conflicts", tileStats.templateLengthStatisticsConflicts_);
            if (!tileStats.templateLengthStatisticsConflicts_)
            {
                xmlWriter.writeElement("Stable", tileStats.templateLengthStatistics_.isStable());
                xmlWriter.writeElement("HighStdDev", tileStats.templateLengthStatistics_.getHighStdDev());
                xmlWriter.writeElement("LowStdDev", tileStats.templateLengthStatistics_.getLowStdDev());
                xmlWriter.writeElement("Max", tileStats.templateLengthStatistics_.getMax());
                xmlWriter.writeElement("Median", tileStats.templateLengthStatistics_.getMedian());
                xmlWriter.writeElement("Min", tileStats.templateLengthStatistics_.getMin());
                xmlWriter.writeElement("Nominal1", alignmentModelName(tileStats.templateLengthStatistics_.getBestModel(0)));
                xmlWriter.writeElement("Nominal2", alignmentModelName(tileStats.templateLengthStatistics_.getBestModel(1)));
                xmlWriter.writeElement("Class1", alignmentClassName(tileStats.templateLengthStatistics_.getBestModel(0)));
                xmlWriter.writeElement("Class2", alignmentClassName(tileStats.templateLengthStatistics_.getBestModel(1)));
            }
        }
    }

    ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, passesFilter ? "Pf" : "Raw")
    {
        // pair stats are stored as read index 0 tiles
        if (dumpReadIndependentStats)
        {
            ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "AlignmentModel")
            {
                for(int model = TemplateLengthStatistics::FFp;
                    TemplateLengthStatistics::InvalidAlignmentModel > model; ++model)
                {
                    xmlWriter.writeElement(alignmentModelName(static_cast<TemplateLengthStatistics::AlignmentModel>(model)), tileStats.alignmentModelCounts_[model]);
                }
                xmlWriter.writeElement("Oversized", tileStats.nominalModelCounts_[TemplateLengthStatistics::Oversized]);
                xmlWriter.writeElement("Undersized", tileStats.nominalModelCounts_[TemplateLengthStatistics::Undersized]);
                xmlWriter.writeElement("Nominal", tileStats.nominalModelCounts_[TemplateLengthStatistics::Nominal]);
                xmlWriter.writeElement("NoMatch", tileStats.nominalModelCounts_[TemplateLengthStatistics::NoMatch]);
            }
            xmlWriter.writeElement("ClusterCount", tileStats.clusterCount_);
            xmlWriter.writeElement("UnanchoredClusterCount", tileStats.unanchoredClusterCount_);
            xmlWriter.writeElement("NmNmClusterCount", tileStats.nmnmClusterCount_);
            xmlWriter.writeElement("RmClusterCount", tileStats.rmClusterCount_);
            xmlWriter.writeElement("QcClusterCount", tileStats.qcClusterCount_);
        }

        ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Read")
        {
            xmlWriter.writeAttribute("number", read);
            ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "AllFragments")
            {
                xmlWriter.writeElement("AlignmentScoreSum", tileStats.alignmentScoreSum_);
                xmlWriter.writeElement("Mismatches", tileStats.mismatches_);
                xmlWriter.writeElement("Count", tileStats.fragmentCount_);
                xmlWriter.writeElement("BasesOutsideIndels", tileStats.basesOutsideIndels_);
                xmlWriter.writeElement("QualityScoreSum", tileStats.qualityScoreSum_);
                xmlWriter.writeElement("Yield", tileStats.yield_);
                xmlWriter.writeElement("YieldQ30", tileStats.yieldQ30_);
                xmlWriter.writeElement("AlignedCount", tileStats.alignedFragmentCount_);
            }

            ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "UniquelyAlignedFragments")
            {
                xmlWriter.writeElement("Mismatches", tileStats.uniquelyAlignedMismatches_);
                xmlWriter.writeElement("Count", tileStats.uniquelyAlignedFragmentCount_);
                xmlWriter.writeElement("Perfect", tileStats.uniquelyAlignedPerfectFragmentCount_);
                xmlWriter.writeElement("BasesOutsideIndels", tileStats.uniquelyAlignedBasesOutsideIndels_);
            }
        }
    }
}


} //namespace matchSelector
} //namespace alignment
} //namespace isaac

