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
 ** \file BamLayout.cpp
 **
 ** Flowcell file locations and such
 **
 ** \author Roman Petrovski
 **/

#include "common/FileSystem.hh"
#include "flowcell/BamLayout.hh"

namespace isaac
{
namespace flowcell
{

template<>
const boost::filesystem::path &Layout::getAttribute<Layout::Bam, BamFilePathAttributeTag>(
    boost::filesystem::path &result) const
{
    ISAAC_ASSERT_MSG(Bam == format_, BamFilePathAttributeTag() << " is only allowed for bam flowcells");

    // boost 1.46 implementation of filesystem::path is coded to instantiate an std::string
    // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());
    pathInternalStringRef = baseCallsPath_.string().c_str();
    return result;
}


} // namespace flowcell
} // namespace isaac
