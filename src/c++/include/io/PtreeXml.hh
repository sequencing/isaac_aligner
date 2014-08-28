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
 ** \file PtreeXml.hh
 **
 ** \brief Simple utility for boost::ptree to XML transformation.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_PTREE_XML_HH
#define iSAAC_IO_PTREE_XML_HH

#include <ostream>
#include <vector>

#include <boost/property_tree/ptree.hpp>

namespace isaac
{
namespace io
{

void index(std::vector<std::string> indexAttributes, boost::property_tree::ptree &tree);

std::ostream &serializeAsXml(std::ostream &os, const boost::property_tree::ptree &tree);

} // namespace io
} // namespace isaac


#endif // iSAAC_IO_PTREE_XML_HH
