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

MatchSelectorStatsXml::MatchSelectorStatsXml(const flowcell::FlowcellLayoutList &flowcellLayoutList)
{
    BOOST_FOREACH(const flowcell::Layout& flowcell, flowcellLayoutList)
    {
        BOOST_FOREACH(const flowcell::ReadMetadata& read, flowcell.getReadMetadataList())
        {
            const std::string readValuePrefix("Stats.<indexed>Flowcell.<flowcell-id>" + flowcell.getFlowcellId() +
                ".<indexed>Read.<number>" + boost::lexical_cast<std::string>(read.getIndex() + 1)
            );
            add(readValuePrefix + ".Length", read.getLength());
        }
    }
}

void MatchSelectorStatsXml::addTile(
    const flowcell::ReadMetadata &read,
    const flowcell::TileMetadata &tile,
    const bool passesFilter,
    const TileStats& tileStats)
{
    const std::string tileValuePrefix("Stats"
                                      ".<indexed>Flowcell.<flowcell-id>" + tile.getFlowcellId()
                                      +".<indexed>Lane.<number>" + boost::lexical_cast<std::string>(tile.getLane())
                                      +".<indexed>Tile.<number>" + boost::lexical_cast<std::string>(tile.getTile())
                                      );

    const std::string tileValuePfPrefix = tileValuePrefix + (passesFilter ? ".Pf" : ".Raw");
    // pair stats are stored as read 0 tiles
    if (!read.getIndex())
    {
        BOOST_STATIC_ASSERT(sizeof(tileStats.alignmentScoreTemplateMismatches_) == sizeof(tileStats.alignmentScoreTemplates_));
        for(unsigned long score = 0; score <= tileStats.maxAlignmentScore_; ++score)
        {
            if (tileStats.alignmentScoreTemplates_[score] || tileStats.alignmentScoreMismatches_[score])
            {
                add<const unsigned long>(tileValuePfPrefix + ".AlignmentScoreDistribution.<indexed>Score.<score>" +
                        boost::lexical_cast<std::string>(score) + ".<xmlattr>.templates",
                        tileStats.alignmentScoreTemplates_[score]);
                add<const unsigned long>(tileValuePfPrefix + ".AlignmentScoreDistribution.<indexed>Score.<score>" +
                        boost::lexical_cast<std::string>(score) + ".<xmlattr>.mismatches",
                        tileStats.alignmentScoreTemplateMismatches_[score]);
            }
        }
    }

    // fragment stats
    const std::string tileReadValuePrefix(
        tileValuePfPrefix +".<indexed>Read.<number>" + boost::lexical_cast<std::string>(read.getIndex() + 1));

    BOOST_STATIC_ASSERT(sizeof(tileStats.alignmentScoreMismatches_) == sizeof(tileStats.alignmentScoreFragments_));
    for(unsigned long score = 0; score <= tileStats.maxAlignmentScore_; ++score)
    {
        if (tileStats.alignmentScoreFragments_[score] || tileStats.alignmentScoreMismatches_[score])
        {
            add<const unsigned long>(tileReadValuePrefix + ".AllFragments.AlignmentScoreDistribution.<indexed>Score.<score>" +
                    boost::lexical_cast<std::string>(score) + ".<xmlattr>.fragments",
                    tileStats.alignmentScoreFragments_[score]);
            add<const unsigned long>(tileReadValuePrefix + ".AllFragments.AlignmentScoreDistribution.<indexed>Score.<score>" +
                    boost::lexical_cast<std::string>(score) + ".<xmlattr>.mismatches",
                    tileStats.alignmentScoreMismatches_[score]);
        }
    }

    for(unsigned long cycle = read.getFirstCycle(); cycle <= read.getLastCycle(); ++cycle)
    {
        /*if (tileStats.cycleBlanks_[cycle] ||
            tileStats.cycleUniquelyAlignedBlanks_[cycle] ||
            tileStats.cycleMismatches_[cycle])*/ //this saves the space but produces xml that the transformations can't deal with well in edge cases
        {
            add<const unsigned long>(tileReadValuePrefix + ".AllFragments.MismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.blanks",
                    tileStats.cycleBlanks_[cycle]);
            add<const unsigned long>(tileReadValuePrefix + ".AllFragments.MismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.mismatches",
                    tileStats.cycleMismatches_[cycle]);
            add<const unsigned long>(tileReadValuePrefix + ".AllFragments.FragmentMismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.one",
                    tileStats.cycle1MismatchFragments_[cycle]);
            add<const unsigned long>(tileReadValuePrefix + ".AllFragments.FragmentMismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.two",
                    tileStats.cycle2MismatchFragments_[cycle]);
            add<const unsigned long>(tileReadValuePrefix + ".AllFragments.FragmentMismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.three",
                    tileStats.cycle3MismatchFragments_[cycle]);
            add<const unsigned long>(tileReadValuePrefix + ".AllFragments.FragmentMismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.four",
                    tileStats.cycle4MismatchFragments_[cycle]);
            add<const unsigned long>(tileReadValuePrefix + ".AllFragments.FragmentMismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.more",
                    tileStats.cycleMoreMismatchFragments_[cycle]);

            add<const unsigned long>(tileReadValuePrefix + ".UniquelyAlignedFragments.MismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.blanks",
                    tileStats.cycleUniquelyAlignedBlanks_[cycle]);
            add<const unsigned long>(tileReadValuePrefix + ".UniquelyAlignedFragments.MismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.mismatches",
                    tileStats.cycleUniquelyAlignedMismatches_[cycle]);
            add<const unsigned long>(tileReadValuePrefix + ".UniquelyAlignedFragments.FragmentMismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.one",
                    tileStats.cycleUniquelyAligned1MismatchFragments_[cycle]);
            add<const unsigned long>(tileReadValuePrefix + ".UniquelyAlignedFragments.FragmentMismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.two",
                    tileStats.cycleUniquelyAligned2MismatchFragments_[cycle]);
            add<const unsigned long>(tileReadValuePrefix + ".UniquelyAlignedFragments.FragmentMismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.three",
                    tileStats.cycleUniquelyAligned3MismatchFragments_[cycle]);
            add<const unsigned long>(tileReadValuePrefix + ".UniquelyAlignedFragments.FragmentMismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.four",
                    tileStats.cycleUniquelyAligned4MismatchFragments_[cycle]);
            add<const unsigned long>(tileReadValuePrefix + ".UniquelyAlignedFragments.FragmentMismatchesByCycle.<indexed>Cycle.<number>" +
                    boost::lexical_cast<std::string>(cycle) + ".<xmlattr>.more",
                    tileStats.cycleUniquelyAlignedMoreMismatchFragments_[cycle]);

        }
    }
    add(tileReadValuePrefix + ".UniquelyAlignedFragments.Count", tileStats.uniquelyAlignedFragmentCount_);
}

void MatchSelectorStatsXml::addTileBarcode(
    const std::string &flowcellId,
    const std::string &projectName, const std::string &sampleName,
    const std::string &barcodeName,
    const flowcell::ReadMetadata &read,
    const flowcell::TileMetadata &tile,
    const bool passesFilter,
    const TileBarcodeStats& tileStats)
{
    const boost::property_tree::path tileValuePrefix("Stats"
                                      "/<indexed>Flowcell/<flowcell-id>" + flowcellId
                                      +"/<indexed>Project/<name>" + projectName
                                      +"/<indexed>Sample/<name>" + sampleName
                                      +"/<indexed>Barcode/<name>" + barcodeName
                                      +"/<indexed>Lane/<number>" + boost::lexical_cast<std::string>(tile.getLane())
                                      +"/<indexed>Tile/<number>" + boost::lexical_cast<std::string>(tile.getTile())
                                      , '/');

    // template length stats stored as read 0 tiles with passesFilter == false
    if (!read.getIndex() && !passesFilter)
    {
        const boost::property_tree::path assumedTemplateLengthPath(tileValuePrefix / "AssumedTemplateLength");
        add(assumedTemplateLengthPath / "Conflicts", tileStats.templateLengthStatisticsConflicts_);
        if (!tileStats.templateLengthStatisticsConflicts_)
        {
            add(assumedTemplateLengthPath / "Stable", tileStats.templateLengthStatistics_.isStable());
            add(assumedTemplateLengthPath / "HighStdDev", tileStats.templateLengthStatistics_.getHighStdDev());
            add(assumedTemplateLengthPath / "LowStdDev", tileStats.templateLengthStatistics_.getLowStdDev());
            add(assumedTemplateLengthPath / "Max", tileStats.templateLengthStatistics_.getMax());
            add(assumedTemplateLengthPath / "Median", tileStats.templateLengthStatistics_.getMedian());
            add(assumedTemplateLengthPath / "Min", tileStats.templateLengthStatistics_.getMin());
            add(assumedTemplateLengthPath / "Nominal1", alignmentModelName(tileStats.templateLengthStatistics_.getBestModel(0)));
            add(assumedTemplateLengthPath / "Nominal2", alignmentModelName(tileStats.templateLengthStatistics_.getBestModel(1)));
            add(assumedTemplateLengthPath / "Class1", alignmentClassName(tileStats.templateLengthStatistics_.getBestModel(0)));
            add(assumedTemplateLengthPath / "Class2", alignmentClassName(tileStats.templateLengthStatistics_.getBestModel(1)));
        }
    }

    const boost::property_tree::path tileValuePfPrefix = tileValuePrefix / (passesFilter ? "Pf" : "Raw");
    // pair stats are stored as read 0 tiles
    if (!read.getIndex())
    {
        for(int model = TemplateLengthStatistics::FFp;
            TemplateLengthStatistics::InvalidAlignmentModel > model; ++model)
        {
            add(tileValuePfPrefix / "AlignmentModel" /
                boost::property_tree::path(alignmentModelName(static_cast<TemplateLengthStatistics::AlignmentModel>(model))),
                tileStats.alignmentModelCounts_[model]);
        }
        add(tileValuePfPrefix / "ClusterCount", tileStats.clusterCount_);
        add(tileValuePfPrefix / "UnanchoredClusterCount", tileStats.unanchoredClusterCount_);
        add(tileValuePfPrefix / "NmNmClusterCount", tileStats.nmnmClusterCount_);
        add(tileValuePfPrefix / "RmClusterCount", tileStats.rmClusterCount_);
        add(tileValuePfPrefix / "QcClusterCount", tileStats.qcClusterCount_);
        add(tileValuePfPrefix / "AlignmentModel" / "Oversized", tileStats.nominalModelCounts_[TemplateLengthStatistics::Oversized]);
        add(tileValuePfPrefix / "AlignmentModel" / "Undersized", tileStats.nominalModelCounts_[TemplateLengthStatistics::Undersized]);
        add(tileValuePfPrefix / "AlignmentModel" / "Nominal", tileStats.nominalModelCounts_[TemplateLengthStatistics::Nominal]);
        add(tileValuePfPrefix / "AlignmentModel" / "NoMatch", tileStats.nominalModelCounts_[TemplateLengthStatistics::NoMatch]);
    }

    // fragment stats
    const boost::property_tree::path tileReadValuePrefix(
        tileValuePfPrefix / "<indexed>Read" / boost::property_tree::path("<number>" + boost::lexical_cast<std::string>(read.getIndex() + 1)));

    add(tileReadValuePrefix / "AllFragments" / "AlignmentScoreSum", tileStats.alignmentScoreSum_);
    add(tileReadValuePrefix / "AllFragments" / "Mismatches", tileStats.mismatches_);
    add(tileReadValuePrefix / "AllFragments" / "Count", tileStats.fragmentCount_);
    add(tileReadValuePrefix / "AllFragments" / "BasesOutsideIndels", tileStats.basesOutsideIndels_);
    add(tileReadValuePrefix / "AllFragments" / "QualityScoreSum", tileStats.qualityScoreSum_);
    add(tileReadValuePrefix / "AllFragments" / "Yield", tileStats.yield_);
    add(tileReadValuePrefix / "AllFragments" / "YieldQ30", tileStats.yieldQ30_);


    add(tileReadValuePrefix / "UniquelyAlignedFragments" / "Mismatches", tileStats.uniquelyAlignedMismatches_);
    add(tileReadValuePrefix / "UniquelyAlignedFragments" / "Count", tileStats.uniquelyAlignedFragmentCount_);
    add(tileReadValuePrefix / "UniquelyAlignedFragments" / "Perfect", tileStats.uniquelyAlignedPerfectFragmentCount_);
    add(tileReadValuePrefix / "UniquelyAlignedFragments" / "BasesOutsideIndels", tileStats.uniquelyAlignedBasesOutsideIndels_);
}

void MatchSelectorStatsXml::addBarcode(
    const flowcell::BarcodeMetadata &barcode)
{
    const std::string barcodeValuePrefix("Stats"
                                      ".<indexed>Flowcell.<flowcell-id>" + barcode.getFlowcellId()
                                      +".<indexed>Barcode.<name>" +  barcode.getName()
                                      +".<indexed>Lane.<number>" + boost::lexical_cast<std::string>(barcode.getLane())
                                      );

    add(barcodeValuePrefix + ".ReferenceName", barcode.getReference());
}

} //namespace matchSelector
} //namespace alignment
} //namespace isaac

