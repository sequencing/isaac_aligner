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

ReadMetadata::ReadMetadata(
    const unsigned number,
    const std::vector<unsigned> &cycleList,
    unsigned index, unsigned
    offset,
    unsigned firstReadCycle)
    : number_(number)
    , cycleList_(cycleList)
    , index_(index)
    , offset_(offset)
    , firstReadCycle_(firstReadCycle)
{
}


ReadMetadata::ReadMetadata(unsigned firstCycle, unsigned lastCycle, unsigned index, unsigned offset)
    : number_(index + 1)
    , cycleList_()
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

unsigned getMaxCycleNumber(const ReadMetadataList &readMetadataList)
{
    using boost::lambda::_1;
    using boost::lambda::_2;
    using boost::lambda::bind;
    return std::accumulate(
        readMetadataList.begin(), readMetadataList.end(), 0,
        bind(std::plus<unsigned>(),
             _1,
             bind<const unsigned>(&flowcell::ReadMetadata::getLastCycle, _2)));
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
