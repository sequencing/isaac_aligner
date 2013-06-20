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
 ** \file FragmentCollector.cpp
 **
 ** \brief See FragmentCollector.hh
 ** 
 ** \author Roman Petrovski
 **/

#include <iostream>
#include <fstream>
#include <cerrno>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "alignment/Cigar.hh"
#include "alignment/matchSelector/FragmentCollector.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "io/Fragment.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

FragmentCollector::~FragmentCollector()
{
}

void FragmentCollector::add(
    const alignment::BamTemplate &bamTemplate,
    const unsigned fragmentIndex,
    const unsigned barcodeIdx)
{
    ISAAC_ASSERT_MSG(2 >= bamTemplate.getFragmentCount(), "Expected paired or single-ended data");

    const alignment::FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(fragmentIndex);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "FragmentCollector::add: " << fragment);
    FragmentBuffer::IndexRecord &recordStart =
        buffer_.initialize(fragment.getCluster().getId(), fragment.getReadIndex());
    recordStart.fStrandPos_ = fragment.getFStrandReferencePosition();
    storeBclAndCigar(fragment, recordStart);

    if (2 == bamTemplate.getFragmentCount())
    {
        const alignment::FragmentMetadata &mate = bamTemplate.getMateFragmentMetadata(fragment);
        if (fragment.isNoMatch())
        {
            ISAAC_ASSERT_MSG(mate.isNoMatch(), "If mate is not a no-match, fragment must be a shadow. fragment: "
                             << fragment << " mate:" << mate);
            recordStart.fragmentHeader() = io::FragmentHeader(bamTemplate, fragment, mate, barcodeIdx, 0);
        }
        else
        {
            const unsigned mateStorageBin = binIndexMap_.getBinIndex(mate.getFStrandReferencePosition());
            recordStart.fragmentHeader() = io::FragmentHeader(bamTemplate, fragment, mate, barcodeIdx,
                                                              mateStorageBin);
        }
    }
    else
    {
        recordStart.fragmentHeader() = io::FragmentHeader(bamTemplate, fragment, barcodeIdx);
    }
}

void FragmentCollector::storeBclAndCigar(
    const alignment::FragmentMetadata & fragment,
    FragmentBuffer::IndexRecord & recordStart)
{
    char *variableData = recordStart.fragmentData();
    // copy the bcl data (reverse-complement the sequence if the fragment is reverse-aligned)
    std::vector<char>::const_iterator bclData = fragment.getBclData();
    if (fragment.isReverse())
    {
        variableData = std::transform(std::reverse_iterator<std::vector<char>::const_iterator>(bclData + fragment.getReadLength()),
                                      std::reverse_iterator<std::vector<char>::const_iterator>(bclData),
                                      variableData, oligo::getReverseBcl);
    }
    else
    {
        variableData = std::copy(bclData, bclData + fragment.getReadLength(), variableData);

    }

    if (fragment.isAligned())
    {
        const alignment::Cigar::const_iterator cigarBegin = fragment.cigarBuffer->begin() + fragment.cigarOffset;
        const alignment::Cigar::const_iterator cigarEnd = cigarBegin + fragment.cigarLength;
        variableData = reinterpret_cast<char*>(std::copy(cigarBegin, cigarEnd, reinterpret_cast<unsigned*>(variableData)));
    }
}


} // namespace matchSelector
} // namespace alignment
} // namespace isaac
