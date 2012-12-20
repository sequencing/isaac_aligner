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
 ** \file DemultiplexingStatsXml.hh
 **
 ** \brief Xml Serialization of Demultiplexing statistics.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_DEMULTIPLEXING_DEMULTIPLEXING_STATS_XML_H
#define ISAAC_DEMULTIPLEXING_DEMULTIPLEXING_STATS_XML_H

#include "io/PtreeXml.hh"
#include "DemultiplexingStats.hh"

namespace isaac
{
namespace demultiplexing
{

class DemultiplexingStatsXml : public boost::property_tree::ptree
{
public:
    DemultiplexingStatsXml();
    void addTileBarcode(
        const std::string &flowcellId,
        const std::string &projectName,
        const std::string &sampleName,
        const std::string &barcodeName,
        const flowcell::TileMetadata &tile,
        const TileBarcodeStats& stat);

    void addFlowcellLane(
        const flowcell::Layout &flowcell,
        const unsigned lane,
        const DemultiplexingStats::LaneBarcodeStats& laneStats);
};

inline std::ostream &operator << (std::ostream &os, const DemultiplexingStatsXml &tree)
{
    return isaac::io::serializeAsXml(os, tree);
}

} //namespace demultiplexing
} //namespace isaac

#endif //ISAAC_DEMULTIPLEXING_DEMULTIPLEXING_STATS_XML_H
