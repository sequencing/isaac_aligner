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
 ** \file RunInfoXml.cpp
 **
 ** Run folder RunInfo.xml helper.
 **
 ** \author Roman Petrovski
 **/

#include <boost/assign.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "common/Debug.hh"
#include "io/PtreeXml.hh"

#include "rta/RunInfoXml.hh"

namespace isaac
{
namespace rta
{

namespace ptree=boost::property_tree;


std::vector<RunInfoXml::ReadInfo> RunInfoXml::getReadInfos() const
{
    const ptree::ptree &readsElements(get_child("RunInfo.Run.Reads"));

    std::vector<RunInfoXml::ReadInfo> ret;
    BOOST_FOREACH (const boost::property_tree::ptree::value_type &readElement, readsElements)
    {
        static const std::string readsElementName("Read");
        if (readElement.first == readsElementName)
        {
            ReadInfo read;
            read.number_ = readElement.second.get<unsigned>("<xmlattr>.Number");
            read.numberOfCycles_ = readElement.second.get<unsigned>("<xmlattr>.NumCycles");
            read.isBarcode_ = 'Y' == readElement.second.get<char>("<xmlattr>.IsIndexedRead");
            ret.push_back(read);
        }
    }

    return ret;
}

std::vector<unsigned> RunInfoXml::getLanes() const
{
    const ptree::ptree &flowcellLayout(get_child("RunInfo.Run.FlowcellLayout"));
    const unsigned laneCount = flowcellLayout.get<unsigned>("<xmlattr>.LaneCount");
    const std::vector<unsigned> ret(boost::make_counting_iterator(1U), boost::make_counting_iterator(laneCount + 1U));
    return ret;
}

std::vector<unsigned> RunInfoXml::getTiles(unsigned lane) const
{
    std::vector<unsigned> ret;
    const ptree::ptree &flowcellLayout(get_child("RunInfo.Run.FlowcellLayout"));
//    const unsigned laneCount = flowcellLayout.get<unsigned>("<xmlattr>.LaneCount");
    const unsigned surfaceCount = flowcellLayout.get<unsigned>("<xmlattr>.SurfaceCount");
    const unsigned swathCount = flowcellLayout.get<unsigned>("<xmlattr>.SwathCount");
    const unsigned tileCount = flowcellLayout.get<unsigned>("<xmlattr>.TileCount");
    const unsigned sectionPerLane = flowcellLayout.get<unsigned>("<xmlattr>.SectionPerLane", 0);
    const unsigned lanePerSection = flowcellLayout.get<unsigned>("<xmlattr>.LanePerSection", 0);


    for (unsigned surface = 1; surface <= surfaceCount; ++surface)
    {
        for (unsigned swath = 1; swath <= swathCount; ++swath)
        {
            for (unsigned tile = 1; tile <= tileCount; ++tile)
            {
                if (sectionPerLane)
                {
                    ISAAC_ASSERT_MSG(lanePerSection, "Unexpected Zero lanePerSection");
                    const unsigned firstSection = 1 + sectionPerLane * ((lane-1)/lanePerSection);
                    for (unsigned section = firstSection; section <= firstSection + sectionPerLane; ++section)
                    {
                        ret.push_back(tile + section * 100 + swath * 1000 + surface * 10000);
                    }
                }
                else
                {
                    ret.push_back(tile + swath * 100 + surface * 1000);
                }
            }
        }
    }

    return ret;
}

std::string RunInfoXml::getFlowcellId() const
{
    const ptree::ptree &flowcell(get_child("RunInfo.Run.Flowcell"));
    return flowcell.get_value<std::string>();
}

std::ostream &operator <<(std::ostream &os, const RunInfoXml &tree)
{
    return io::serializeAsXml(os, tree);
}

std::istream &operator >> (std::istream &is, RunInfoXml &tree)
{
    boost::property_tree::read_xml<boost::property_tree::ptree>(is, tree);
    return is;
}

} // namespace rta
} // namespace isaac

