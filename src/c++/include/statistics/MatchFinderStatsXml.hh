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
