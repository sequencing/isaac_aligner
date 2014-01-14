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
 ** \file RunInfoXml.hh
 **
 ** Run folder RunInfo.xml helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_RTA_RUN_INFO_XML_HH
#define iSAAC_RTA_RUN_INFO_XML_HH

#include <vector>
#include <utility>
#include <string>
#include <boost/property_tree/ptree.hpp>

namespace isaac
{
namespace rta
{

class RunInfoXml : public boost::property_tree::ptree
{
public:
    struct ReadInfo
    {
        unsigned number_;
        unsigned numberOfCycles_;
        bool isBarcode_;
    };
    std::vector<ReadInfo> getReadInfos() const;
    std::vector<unsigned> getLanes() const;
    std::vector<unsigned> getTiles(unsigned lane) const;
    std::string getFlowcellId() const;
};

std::ostream &operator << (std::ostream &os, const RunInfoXml &tree);
std::istream &operator >> (std::istream &is, RunInfoXml &tree);

} // namespace rta
} // namespace isaac

#endif // #ifndef iSAAC_RTA_RUN_INFO_XML_HH
