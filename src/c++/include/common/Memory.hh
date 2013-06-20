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

static const unsigned long PAGE_SIZE = boost::interprocess::mapped_region::get_page_size();

inline unsigned long pageRoundUp(unsigned long size)
{
    return (size + PAGE_SIZE - 1) & (~(PAGE_SIZE - 1));
}


} //namespace common
} //namespace isaac

#endif // #ifndef iSAAC_COMMON_MEMORY_HPP
