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
 ** \file BinIndexMap.hh
 **
 ** \brief Fragment buffer flushing and output file management.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BIN_INDEX_MAP_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_BIN_INDEX_MAP_HH

#include <boost/foreach.hpp>

#include "alignment/MatchDistribution.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

class BinIndexMap;
std::ostream& operator << (std::ostream& os, const BinIndexMap &binIndexMap);
/**
 ** \brief associate a genomic location to a bin index
 **
 ** The map has exactly the same geometry as the original MatchDistribution
 ** (outer vector for the contigs, inner vector for the small bins used to
 ** compute the distribution over each contig). The values stored provide
 ** direct access to the output bin to use for a given reference position.
 **/
class BinIndexMap : public std::vector<std::vector<unsigned> >
{
    /// the binSize from the MatchDistribution
    const unsigned binSize_;
public:
    BinIndexMap(const MatchDistribution &matchDistribution,
                const unsigned long outputBinSize)
        : binSize_(matchDistribution.getBinSize())
    {
        reserve(matchDistribution.size());
        size_t currentBinIndex = 0;
        // first push the contig and bin for unaligned clusters
        push_back(std::vector<unsigned>());
        back().reserve(1);
        back().push_back(currentBinIndex++);

        // now put in all the contig bins
        BOOST_FOREACH(const std::vector<unsigned> &contigDistribution, matchDistribution)
        {
            push_back(std::vector<unsigned>());
            back().reserve(contigDistribution.size());
            unsigned long currentBinSize = 0;
            BOOST_FOREACH(unsigned count, contigDistribution)
            {
                // the indexes produced are expected to be contiguous.
                // avoid advancing bin index if it has not been pushed even
                // if count already exceeds outputBinSize
                if (currentBinSize && currentBinSize + count > outputBinSize)
                {
                    ++currentBinIndex;
                    currentBinSize = 0;
                }
                currentBinSize += count;
                back().push_back(currentBinIndex);
            }
            // bins do not spread across contigs
            ++currentBinIndex;
        }
    }

    /**
     ** \brief convert a reference position on a contig into a bin index that
     ** can be used to identify either the file path or the stream associated
     ** to the ReferencePosition.
     **/
    size_t getBinIndex(const reference::ReferencePosition &referencePosition) const
    {
        const unsigned long contigId = referencePosition.getContigId();
        const std::vector<unsigned> &binIndexList = at(contigId + 1);
        const unsigned long position = referencePosition.getPosition();
        const unsigned long index = position / binSize_;
        assert(binIndexList.size() > index);
        return binIndexList[index];
    }

    /**
     * \return The first reference position that can be found in the bin
     */
    reference::ReferencePosition getBinFirstPos(size_t bin) const
    {
        std::vector<unsigned>::const_reference (std::vector<unsigned>::*b)() const = &std::vector<unsigned>::back;

        const_iterator binContigIterator = lower_bound(begin(), end(), bin, boost::bind(b, _1) < _2);

        ISAAC_ASSERT_MSG(end() != binContigIterator, "Bin number has to be one of those we have a contig for");
        std::vector<unsigned>::const_iterator contigBinIterator =
            lower_bound(binContigIterator->begin(), binContigIterator->end(), bin);

        ISAAC_ASSERT_MSG(binContigIterator->end() != contigBinIterator, "Bin number must be present in the contig bins");

        const unsigned long contigId = binContigIterator - begin() - 1;
        const unsigned long position = binSize_ * (contigBinIterator - binContigIterator->begin());

        return reference::ReferencePosition(contigId, position);
    }

    /**
     * \return The first reference position that belongs to the subsequent bin. NOTE: for last bin in the contig
     *         there is not guarantee that no alignments will exist at this position and beyond. However, the amount
     *         of data aligning there should be considered minor and belonging to the last bin.
     */
    reference::ReferencePosition getBinFirstInvalidPos(size_t bin) const
    {
        std::vector<unsigned>::const_reference (std::vector<unsigned>::*b)() const = &std::vector<unsigned>::back;

        const_iterator binContigIterator = lower_bound(begin(), end(), bin, boost::bind(b, _1) < _2);

        ISAAC_ASSERT_MSG(end() != binContigIterator, "Bin number has to be one of those we have a contig for");

        std::vector<unsigned>::const_iterator contigNextBinIterator =
            upper_bound(binContigIterator->begin(), binContigIterator->end(), bin);

        const unsigned long contigId = binContigIterator - begin() - 1;
        const unsigned long position = binSize_ * (contigNextBinIterator - binContigIterator->begin());

        return reference::ReferencePosition(contigId, position);
    }

    unsigned getBinCount() const
    {
        return back().back() + 1;
    }

};


inline std::ostream& operator << (std::ostream& os, const BinIndexMap &binIndexMap)
{
    std::string message;
    BOOST_FOREACH(const std::vector<unsigned> &contigIndexList, binIndexMap)
    {
        message += contigIndexList.empty() ?
            std::string("Empty bin list") :
            (boost::format("%d bin indexes from %d to %d:") % contigIndexList.size() % contigIndexList.front() % contigIndexList.back()).str();
/*
        BOOST_FOREACH(unsigned index, contigIndexList)
        {
            message += (boost::format(" %d") % index).str();
        }
*/
        message += "\n";
    }
    os << message;
    return os;
}

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BIN_INDEX_MAP_HH
