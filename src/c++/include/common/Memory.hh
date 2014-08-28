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
 ** \file Memory.hh
 **
 ** memory management helper utilities.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_COMMON_MEMORY_HPP
#define iSAAC_COMMON_MEMORY_HPP

#include <boost/interprocess/mapped_region.hpp>

#include "common/Debug.hh"

namespace isaac
{
namespace common
{
// PAGE_SIZE is taken by some macro in cygwin...
static const unsigned long ISAAC_PAGE_SIZE = boost::interprocess::mapped_region::get_page_size();

inline unsigned long pageRoundUp(unsigned long size)
{
    return (size + ISAAC_PAGE_SIZE - 1) & (~(ISAAC_PAGE_SIZE - 1));
}


} //namespace common
} //namespace isaac

#endif // #ifndef iSAAC_COMMON_MEMORY_HPP
