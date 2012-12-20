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
    ISAAC_THREAD_CERR_DEV_TRACE((boost::format("FragmentCollector::add: %s") % fragment).str());
    FragmentBuffer::IndexRecord &recordStart =
        buffer_.initialize(fragment.getCluster().getId(), fragment.getReadIndex());
    recordStart.fStrandPos_ = fragment.getFStrandReferencePosition();
    storeBclAndCigar(fragment, recordStart);

    if (2 == bamTemplate.getFragmentCount())
    {
        const alignment::FragmentMetadata &mate = bamTemplate.getMateFragmentMetadata(fragment);
        if (mate.isNoMatch() && fragment.isNoMatch())
        {
            recordStart.nmIndex() = io::NmFragmentIndex();
            recordStart.fragmentHeader() = io::FragmentHeader(bamTemplate, fragment, mate, barcodeIdx,
                                                              true);
        }
        else
        {
            const unsigned mateStorageBin = binIndexMap_.getBinIndex(mate.getFStrandReferencePosition());

            if (!(fragment.isReverse() || !fragment.isAligned()))
            {
                recordStart.fIndex() = io::FStrandFragmentIndex(
                    fragment.getFStrandReferencePosition(),
                    io::FragmentIndexMate(mate, mateStorageBin),
                    bamTemplate);
            }
            else
            {
                recordStart.rsIndex() = io::RStrandOrShadowFragmentIndex(
                    (fragment.isAligned() ? fragment.getFStrandReferencePosition()
                        // shadows are stored at the position of their singletons
                        : mate.getFStrandReferencePosition()),
                          io::FragmentIndexAnchor(fragment),
                          io::FragmentIndexMate(mate, mateStorageBin),
                          bamTemplate);
            }
            const unsigned storageBin = binIndexMap_.getBinIndex(fragment.getFStrandReferencePosition());
            recordStart.fragmentHeader() = io::FragmentHeader(bamTemplate, fragment, mate, barcodeIdx,
                                                              storageBin == mateStorageBin);
        }
    }
    else
    {
        if (!fragment.isNoMatch())
        {
            recordStart.seIndex() = io::SeFragmentIndex(fragment.getFStrandReferencePosition());
        }
        else
        {
            recordStart.nmIndex() = io::NmFragmentIndex();
        }

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
    const alignment::Cigar::const_iterator cigarBegin = fragment.cigarBuffer->begin() + fragment.cigarOffset;
    const alignment::Cigar::const_iterator cigarEnd = cigarBegin + fragment.cigarLength;
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
    variableData = reinterpret_cast<char*>(std::copy(cigarBegin, cigarEnd, reinterpret_cast<unsigned*>(variableData)));
}


} // namespace matchSelector
} // namespace alignment
} // namespace isaac
