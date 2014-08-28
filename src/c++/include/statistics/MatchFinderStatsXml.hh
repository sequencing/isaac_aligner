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
 ** \file MatchFinderStatsXml.hh
 **
 ** \brief Xml Serialization of MatchFinder statistics.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_STATISTICS_MATCH_FINDER_STATS_XML_H
#define ISAAC_STATISTICS_MATCH_FINDER_STATS_XML_H

#include <ostream>

#include "io/PtreeXml.hh"
#include "statistics/MatchFinderTileStats.hh"

namespace isaac
{
namespace statistics{

class MatchFinderStatsXml : public boost::property_tree::ptree
{
public:
    void addTile(const std::string &sampleName, const std::string &barcode,
                 const unsigned seedOffset,
                 const unsigned lane, const unsigned read,
                 const unsigned tile, const MatchFinderTileStats& tileStats);
};

inline std::ostream &operator << (std::ostream &os, const MatchFinderStatsXml &tree)
{
    return isaac::io::serializeAsXml(os, tree);
}


} //namespace statistics
} //namespace isaac

#endif //ISAAC_STATISTICS_MATCH_FINDER_STATS_XML_H
