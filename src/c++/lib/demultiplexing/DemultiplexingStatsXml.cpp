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
 ** \file DemultiplexingStatsXml.cpp
 **
 ** \brief Xml Serialization of Demultiplexing statistics.
 **
 ** \author Roman Petrovski
 **/

#include <boost/lexical_cast.hpp>

#include "demultiplexing/DemultiplexingStatsXml.hh"

namespace isaac
{
namespace demultiplexing
{

DemultiplexingStatsXml::DemultiplexingStatsXml()
{
}

void DemultiplexingStatsXml::addTileBarcode(
    const std::string &flowcellId,
    const std::string &projectName, const std::string &sampleName,
    const std::string &barcodeName,
    const flowcell::TileMetadata &tile,
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

    add(tileValuePrefix / "BarcodeCount", tileStats.barcodeCount_);
    add(tileValuePrefix / "PerfectBarcodeCount", tileStats.perfectBarcodeCount_);
    add(tileValuePrefix / "OneMismatchBarcodeCount", tileStats.oneMismatchBarcodeCount_);
}

void DemultiplexingStatsXml::addFlowcellLane(
    const flowcell::Layout &flowcell,
    const unsigned lane,
    const DemultiplexingStats::LaneBarcodeStats& laneStats)
{
    const boost::property_tree::path laneValuePrefix("Stats"
                                      "/<indexed>Flowcell/<flowcell-id>" + flowcell.getFlowcellId()
                                      +"/<indexed>Lane/<number>" + boost::lexical_cast<std::string>(lane)
                                      +"/TopUnknownBarcodes", '/');

    BOOST_FOREACH(const DemultiplexingStats::UnknownBarcodeHits::value_type &unknownBarcode,
                  laneStats.topUnknownBarcodes_)
    {
        add(laneValuePrefix / ("<indexed>Barcode/<sequence>" + bases(unknownBarcode.first, flowcell.getBarcodeLength())).c_str()
            / "<xmlattr>" / "count", unknownBarcode.second);
    }
}

} //namespace demultiplexing
} //namespace isaac

