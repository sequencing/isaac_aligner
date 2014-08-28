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
 ** \file BoostArchiveHelpers.hh
 **
 ** \brief helpers for serializing various basic types
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_COMMON_BOOST_ARCHIVE_HELPERS_H
#define ISAAC_COMMON_BOOST_ARCHIVE_HELPERS_H

#include <boost/serialization/split_free.hpp>

#include "alignment/MatchTally.hh"
#include "workflow/AlignWorkflow.hh"

BOOST_SERIALIZATION_SPLIT_FREE(boost::filesystem::path)

/**
 * \brief serialization implementation types that don't require private member access
 * can be placed in the boost::serialization namespace
 */
namespace boost {
namespace serialization {

template <class Archive>
void save(Archive &ar, const boost::filesystem::path &p, const unsigned int version)
{
    ar << boost::serialization::make_nvp("path", p.string());
}

template <class Archive>
void load(Archive &ar, boost::filesystem::path &p, const unsigned int version)
{
    std::string tmp;
    ar >> boost::serialization::make_nvp("path", tmp);
    p = tmp;
}

template <class Archive>
void serialize(Archive &ar, std::pair<std::string, std::string> &ps, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(ps.first);
    ar & BOOST_SERIALIZATION_NVP(ps.second);
}

} //namespace serialization
} //namespace boost

#endif //ISAAC_COMMON_BOOST_ARCHIVE_HELPERS_H
