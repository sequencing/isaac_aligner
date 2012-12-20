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
 ** \file MatchFinderStatsXmlWriter.cpp
 **
 ** \brief Xml Serialization of MatchFinder statistics.
 **
 ** \author Roman Petrovski
 **/

#include <boost/lexical_cast.hpp>

#include "statistics/MatchFinderStatsXml.hh"


namespace isaac
{
namespace statistics
{

void MatchFinderStatsXml::addTile(
    const std::string &sampleName, const std::string &barcode,
    const unsigned read, const unsigned seedOffset,
    const unsigned lane, const unsigned tile,
    const MatchFinderTileStats& tileStats)
{
    const std::string tileValuePrefix("Stats"
                                      ".<indexed>Sample.<name>" + sampleName
                                      +".<indexed>Barcode.<name>" + barcode
                                      +".<indexed>Read.<number>" + boost::lexical_cast<std::string>(read)
                                      +".<indexed>Seed.<offset>" + boost::lexical_cast<std::string>(seedOffset)
                                      +".<indexed>Lane.<number>" + boost::lexical_cast<std::string>(lane)
                                      +".<indexed>Tile.<number>" + boost::lexical_cast<std::string>(tile)
                                      );

    add(tileValuePrefix + ".NoMatchSeeds", tileStats.noMatchSeeds_);
    add(tileValuePrefix + ".RepeatMatchSeeds", tileStats.repeatMatchSeeds_);
    add(tileValuePrefix + ".RepeatMatches", tileStats.repeatMatches_);
    add(tileValuePrefix + ".UniqueMatchSeeds", tileStats.uniqueMatchSeeds_);
    add(tileValuePrefix + ".TooManyRepeatsSeeds", tileStats.tooManyRepeatsSeeds_);
}


} //namespace statistics
} //namespace isaac

