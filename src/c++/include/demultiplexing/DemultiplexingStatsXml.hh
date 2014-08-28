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
    void addLaneBarcode(
        const std::string &flowcellId,
        const std::string &projectName,
        const std::string &sampleName,
        const std::string &barcodeName,
        const unsigned lane,
        const LaneBarcodeStats& stat);

    void addFlowcellLane(
        const flowcell::Layout &flowcell,
        const unsigned lane,
        const LaneBarcodeStats& laneStats);
};

inline std::ostream &operator << (std::ostream &os, const DemultiplexingStatsXml &tree)
{
    return isaac::io::serializeAsXml(os, tree);
}

} //namespace demultiplexing
} //namespace isaac

#endif //ISAAC_DEMULTIPLEXING_DEMULTIPLEXING_STATS_XML_H
