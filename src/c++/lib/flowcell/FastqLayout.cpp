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
 ** \file FastqLayout.cpp
 **
 ** Flowcell file locations and such
 **
 ** \author Roman Petrovski
 **/

#include "common/FileSystem.hh"
#include "flowcell/FastqLayout.hh"

namespace isaac
{
namespace flowcell
{

namespace fastq
{

void getFastqFilePath(
    const boost::filesystem::path &baseCallsPath,
    const unsigned lane,
    const unsigned read,
    const bool compressed,
    boost::filesystem::path &result)
{
    ISAAC_ASSERT_MSG(read <= fastq::READ_NUMBER_MAX, "Read number " << read << " must not exceed " << fastq::READ_NUMBER_MAX);

    // Warning: all this mad code below is to avoid memory allocations during path formatting.
    // the result is expected to be pre-sized, else allocations will occur as usual.
    char fastqFileName[100];
    const int fastqFileNameLength = snprintf(fastqFileName, sizeof(fastqFileName),
                                          (compressed ?
                                              "%clane%d_read%d.fastq.gz" : "%clane%d_read%d.fastq"),
                                              common::getDirectorySeparatorChar(), lane, read);

    // boost 1.46 implementation of filesystem::path is coded to instantiate an std::string
    // when doing append. Therefore have to jump through hoops to prevent memory allocations from happening
    std::string & pathInternalStringRef = const_cast<std::string&>(result.string());
    pathInternalStringRef = baseCallsPath.string().c_str();
    pathInternalStringRef.append(fastqFileName, fastqFileName + fastqFileNameLength);
}

} //namespace fastq

template<>
const boost::filesystem::path &Layout::getLaneReadAttribute<Layout::Fastq, FastqFilePathAttributeTag>(
    const unsigned lane, const unsigned read, boost::filesystem::path &result) const
{
    ISAAC_ASSERT_MSG(Fastq == format_, FastqFilePathAttributeTag() << " is only allowed for bam flowcells");
    ISAAC_ASSERT_MSG(lane <= laneNumberMax_, "Lane number " << lane << " must not exceed " << laneNumberMax_);

    const FastqFlowcellData &data = boost::get<FastqFlowcellData>(formatSpecificData_);

    fastq::getFastqFilePath(getBaseCallsPath(), lane, read, data.compressed_, result);
    return result;
}


} // namespace flowcell
} // namespace isaac
