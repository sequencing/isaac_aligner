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
 ** \file MatchSelectorStatsXml.hh
 **
 ** \brief Xml Serialization of MatchSelector statistics.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_ALIGNMENT_MATCH_SELECTOR_STATS_XML_H
#define ISAAC_ALIGNMENT_MATCH_SELECTOR_STATS_XML_H

#include "io/PtreeXml.hh"
#include "MatchSelectorStats.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

class MatchSelectorStatsXml : public boost::property_tree::ptree
{
public:
    MatchSelectorStatsXml(const flowcell::FlowcellLayoutList &flowcellLayoutList);
    void addTile(
        const flowcell::ReadMetadata &read,
        const flowcell::TileMetadata &tile,
        const bool passesFilter,
        const TileStats& tileStats);
    void addTileBarcode(
        const std::string &flowcellId,
        const std::string &projectName,
        const std::string &sampleName,
        const std::string &barcodeName,
        const flowcell::ReadMetadata &read,
        const flowcell::TileMetadata &tile,
        const bool passesFilter,
        const TileBarcodeStats& tileStats);
    void addBarcode(
        const flowcell::BarcodeMetadata &barcode);
};

inline std::ostream &operator << (std::ostream &os, const MatchSelectorStatsXml &tree)
{
    return isaac::io::serializeAsXml(os, tree);
}

} //namespace matchSelector
} //namespace alignment
} //namespace isaac

#endif //ISAAC_ALIGNMENT_MATCH_SELECTOR_STATS_XML_H
