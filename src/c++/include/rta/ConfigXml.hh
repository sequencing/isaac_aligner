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
 ** \file ConfigXml.hh
 **
 ** BaseCalls config.xml helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_RTA_CONFIG_XML_HH
#define iSAAC_RTA_CONFIG_XML_HH

#include <vector>
#include <utility>
#include <string>
#include <boost/property_tree/ptree.hpp>

namespace isaac
{
namespace rta
{

class ConfigXml : public boost::property_tree::ptree
{
public:
    struct RunParametersRead
    {
        unsigned index_;
        unsigned firstCycle_;
        unsigned lastCycle_;
    };
    std::vector<RunParametersRead> getRunParametersReads() const;
    std::vector<unsigned> getLanes() const;
    std::vector<unsigned> getTiles(unsigned lane) const;
    std::pair<std::string, std::string> getSoftwareVersion() const;
    std::string getFlowcellId() const;
};

std::ostream &operator << (std::ostream &os, const ConfigXml &tree);
std::istream &operator >> (std::istream &is, ConfigXml &indexedTree);

} // namespace basecalls
} // namespace isaac

#endif // #ifndef iSAAC_RTA_CONFIG_XML_HH
