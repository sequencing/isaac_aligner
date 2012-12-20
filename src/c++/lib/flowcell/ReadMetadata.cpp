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
 ** \file readMetadata.cpp
 **
 ** Packaging of the metadata associated to a read.
 **
 ** \author Come Raczy
 **/

#include <numeric>
#include <boost/foreach.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "flowcell/ReadMetadata.hh"

namespace isaac
{
namespace flowcell
{

ReadMetadata::ReadMetadata(const std::vector<unsigned> &cycleList, unsigned index, unsigned offset, unsigned firstReadCycle)
    : cycleList_(cycleList)
    , index_(index)
    , offset_(offset)
    , firstReadCycle_(firstReadCycle)
{
}


ReadMetadata::ReadMetadata(unsigned firstCycle, unsigned lastCycle, unsigned index, unsigned offset)
    : cycleList_()
    , index_(index)
    , offset_(offset)
    , firstReadCycle_(firstCycle)
{
    for (unsigned cycle = firstCycle; lastCycle >= cycle; ++cycle)
    {
        cycleList_.push_back(cycle);
    }
}

bool ReadMetadata::operator==(const ReadMetadata &rhs) const
{
    return cycleList_ == rhs.cycleList_ &&
        index_ == index_ &&
        offset_ == offset_;
}

unsigned getTotalReadLength(const ReadMetadataList &readMetadataList)
{
    using boost::lambda::_1;
    using boost::lambda::_2;
    using boost::lambda::bind;
    return std::accumulate(
        readMetadataList.begin(), readMetadataList.end(), 0,
        bind(std::plus<unsigned>(),
             _1,
             bind<const unsigned>(&flowcell::ReadMetadata::getLength, _2)));
}

std::vector<unsigned> getAllCycleNumbers(const ReadMetadataList &readMetadataList)
{
    std::vector<unsigned> cycleNumbers;
    BOOST_FOREACH(const flowcell::ReadMetadata &readMetadata, readMetadataList)
    {
        cycleNumbers.insert(cycleNumbers.end(),
                            boost::make_counting_iterator(readMetadata.getFirstCycle()),
                            boost::make_counting_iterator(readMetadata.getLastCycle() + 1));
    }
    return cycleNumbers;
}


} // namespace flowcell
} // namespace isaac
