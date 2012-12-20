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
 ** \file ConfigXml.cpp
 **
 ** BaseCalls config.xml helper.
 **
 ** \author Roman Petrovski
 **/

#include <boost/assign.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "io/PtreeXml.hh"

#include "basecalls/ConfigXml.hh"

namespace isaac
{
namespace basecalls
{

namespace ptree=boost::property_tree;

std::vector<ConfigXml::RunParametersRead> ConfigXml::getRunParametersReads() const
{
    static const std::string lanesKey("BaseCallAnalysis.Run.RunParameters");
    std::vector<RunParametersRead> ret;
    const ptree::ptree &readsElements(get_child(lanesKey));
    BOOST_FOREACH(const ptree::ptree::value_type &potentiallyReadsNode, readsElements)
    {
        static const std::string index("Reads");
        if(index == potentiallyReadsNode.first)
        {
            RunParametersRead read;
            read.index_ = potentiallyReadsNode.second.get<unsigned>("<xmlattr>.Index");
            read.firstCycle_ = potentiallyReadsNode.second.get<unsigned>("FirstCycle");
            read.lastCycle_ = potentiallyReadsNode.second.get<unsigned>("LastCycle");
            ret.push_back(read);
        }
    }
    return ret;
}


std::vector<unsigned> ConfigXml::getLanes() const
{
    static const std::string lanesKey("BaseCallAnalysis.Run.TileSelection.<indexed>Lane");
    std::vector<unsigned> ret;
    const ptree::ptree &laneElements(get_child(lanesKey));
    BOOST_FOREACH(const ptree::ptree::value_type &laneIndexNode, laneElements)
    {
        static const std::string index("<Index>");
        BOOST_ASSERT(boost::starts_with(laneIndexNode.first, index));
        ret.push_back(boost::lexical_cast<unsigned>(laneIndexNode.first.substr(index.length())));
    }
    std::sort(ret.begin(), ret.end());
    return ret;
}

std::vector<unsigned> ConfigXml::getTiles(unsigned lane) const
{
    const std::string laneKey("BaseCallAnalysis.Run.TileSelection.<indexed>Lane.<Index>" +
                              boost::lexical_cast<std::string>(lane));
    std::vector<unsigned> ret;
    const ptree::ptree &laneElements(get_child(laneKey));
    BOOST_FOREACH(const ptree::ptree::value_type &potentiallyTileNode, laneElements)
    {
        if ("Tile" == potentiallyTileNode.first)
        {
            ret.push_back(potentiallyTileNode.second.get_value<unsigned>());
        }
    }
    std::sort(ret.begin(), ret.end());
    return ret;
}

std::pair<std::string, std::string> ConfigXml::getSoftwareVersion() const
{
    return std::pair<std::string, std::string>(
        get<std::string>("BaseCallAnalysis.Run.Software.<xmlattr>.Name"),
        get<std::string>("BaseCallAnalysis.Run.Software.<xmlattr>.Version"));        
}

std::string ConfigXml::getFlowcellId() const
{
    const ptree::ptree &runParameters = get_child("BaseCallAnalysis.Run.RunParameters");
    const std::string path = "RunFlowcellId";
    if (runParameters.not_found() != runParameters.find(path))
    {
        return runParameters.get<std::string>(path);
    }
//    const std::string path = "BaseCallAnalysis.Run.RunParameters.RunFlowcellId";
//    if (not_found() !=  find(path))
//    {
//        return get<std::string>(path);
//    }
    return "";
}

std::ostream &operator <<(std::ostream &os, const ConfigXml &tree)
{
    return io::serializeAsXml(os, tree);
}

std::istream &operator >> (std::istream &is, ConfigXml &indexedTree)
{
    boost::property_tree::read_xml<boost::property_tree::ptree>(is, indexedTree);
    static const std::vector<std::string > indexAttrs = boost::assign::list_of
            (std::string("BaseCallAnalysis.Run.TileSelection.Lane.Index"))
//            (std::string("BaseCallAnalysis.Run.RunParameters.Reads.Index"))
            ;
    io::index(indexAttrs, indexedTree);
    return is;
}

} // namespace basecalls
} // namespace isaac

