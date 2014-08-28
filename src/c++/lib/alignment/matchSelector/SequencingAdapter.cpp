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
 ** \file SequencingAdapter.cpp
 **
 ** \brief See SequencingAdapter.hh
 ** 
 ** \author Roman Petrovski
 **/
#include <boost/format.hpp>

#include "alignment/matchSelector/SequencingAdapter.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

SequencingAdapter::SequencingAdapter(const flowcell::SequencingAdapterMetadata &adapterMetadata) :
    adapterMetadata_(adapterMetadata),
    kmerPositions_(oligo::getMaxKmer<unsigned>(adapterMatchBasesMin_) + 1, char(UNINITIALIZED_POSITION))
{

    ISAAC_ASSERT_MSG(adapterMetadata_.getSequence().size() < unsigned(std::numeric_limits<char>::max()), "Adapter sequence is too long");
    ISAAC_ASSERT_MSG((adapterMetadata_.isUnbounded() ||
                     adapterMetadata_.getSequence().size() <= adapterMetadata_.getClipLength()),
                     "Clip length cannot be shorter than the adapter sequence");

    oligo::KmerGenerator<unsigned short, std::string::const_iterator> kmerGenerator(
        adapterMetadata_.getSequence().begin(), adapterMetadata_.getSequence().end(), adapterMatchBasesMin_);
    std::string::const_iterator position = adapterMetadata_.getSequence().begin();
    unsigned short kmer = 0;
    while(kmerGenerator.next(kmer, position))
    {
        char &pos = kmerPositions_.at(kmer);
        if (UNINITIALIZED_POSITION == pos)
        {
            pos = std::distance(adapterMetadata_.getSequence().begin(), position);
        }
        else if (NON_UNIQUE_KMER_POSITION != pos)
        {
            pos = NON_UNIQUE_KMER_POSITION;
        }
    }
}

const std::pair<std::vector<char>::const_iterator, std::vector<char>::const_iterator>
SequencingAdapter::getMatchRange(
    const std::vector<char>::const_iterator sequenceBegin,
    const std::vector<char>::const_iterator sequenceEnd,
    const std::vector<char>::const_iterator mismatchBase) const
{
    unsigned short kmer = 0;
    if (oligo::generateKmer(adapterMatchBasesMin_, kmer, mismatchBase, sequenceEnd))
    {
        char pos = kmerPositions_[kmer];
        if (isGoodPosition(pos))
        {
            const unsigned mismatchBaseOffset = std::distance(sequenceBegin, mismatchBase);
            const unsigned adapterBasesBeforeSequence = mismatchBaseOffset < unsigned(pos) ? pos - mismatchBaseOffset : 0;
            if (!adapterBasesBeforeSequence || !adapterMetadata_.isUnbounded())
            {
                const std::vector<char>::const_iterator testBase = mismatchBase - (pos - adapterBasesBeforeSequence);
                const unsigned testSequenceLength = std::distance(testBase, sequenceEnd);
                const unsigned adapterSequenceSize = adapterMetadata_.getSequence().size();
                const unsigned leftClippedAdapaterLength = adapterSequenceSize - adapterBasesBeforeSequence;
                const unsigned overlapLength = std::min<unsigned>(testSequenceLength, leftClippedAdapaterLength);
                if (overlapLength < leftClippedAdapaterLength && adapterMetadata_.isUnbounded() && adapterMetadata_.isReverse())
                {
                    ISAAC_THREAD_CERR_DEV_TRACE("SequencingAdapter::checkSequence: unbounded adapter begins after the reverse sequence. Ignoring: " <<
                                                std::string(mismatchBase, sequenceEnd) <<
                                                " adapterBasesBeforeSequence:" << adapterBasesBeforeSequence <<
                                                " overlapLength:" << leftClippedAdapaterLength <<
                                                " &*testBase:" << std::string(&*testBase, &*testBase + overlapLength));
                }
                else
                {
                    if (overlapLength >= adapterMatchBasesMin_ &&
                        !adapterMetadata_.getSequence().compare(adapterBasesBeforeSequence, overlapLength, &*testBase, overlapLength))
                    {
                        ISAAC_THREAD_CERR_DEV_TRACE("SequencingAdapter::checkSequence: found: " << std::string(mismatchBase, sequenceEnd) <<
                                                    " unbounded: " << adapterMetadata_.isUnbounded() <<
                                                    " adapterBasesBeforeSequence:" << adapterBasesBeforeSequence <<
                                                    " compareFullLength:" << overlapLength <<
                                                    " &*testBase:" << std::string(&*testBase, &*testBase + overlapLength));
                        if (adapterMetadata_.isReverse())
                        {
                            return adapterMetadata_.isUnbounded() ?
                                std::make_pair(sequenceBegin, testBase + overlapLength) :
                                std::make_pair(testBase - std::min<unsigned>(std::distance(sequenceBegin, testBase),
                                                                   adapterMetadata_.getClipLength() - adapterSequenceSize),
                                               testBase + overlapLength);
                        }
                        else
                        {
                            return adapterMetadata_.isUnbounded() ?
                                std::make_pair(testBase, sequenceEnd) :
                                std::make_pair(testBase, testBase + std::min(overlapLength, adapterMetadata_.getClipLength()));
                        }
                    }
                    else
                    {
        //                ISAAC_THREAD_CERR_DEV_TRACE((boost::format("SequencingAdapter::checkSequence: does not compare: %s") %
        //                    std::string(mismatchBase, sequenceEnd)).str());
                    }
                }
            }
            else
            {
                ISAAC_THREAD_CERR_DEV_TRACE("SequencingAdapter::checkSequence: unbounded adapter begins before the forward sequence. Ignoring: " <<
                                            std::string(mismatchBase, sequenceEnd) <<
                                            " adapterBasesBeforeSequence:" << adapterBasesBeforeSequence);
            }

        }
        else
        {
//            ISAAC_THREAD_CERR_DEV_TRACE((boost::format("SequencingAdapter::checkSequence: not found: %s") %
//                std::string(mismatchBase, sequenceEnd)).str());
        }
    }
    else
    {
//        ISAAC_THREAD_CERR_DEV_TRACE((boost::format("SequencingAdapter::checkSequence: too short: %s") %
//            std::string(mismatchBase, sequenceEnd)).str());
    }
    return std::make_pair(mismatchBase, mismatchBase);
}

} // namespace matchSelector
} // namespace alignment
} // namespace isaac
