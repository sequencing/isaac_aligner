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
    const unsigned distributionBinSize_;
public:
    BinIndexMap(const MatchDistribution &matchDistribution,
                const unsigned long outputBinSize,
                const bool skipEmptyBins)
        : distributionBinSize_(matchDistribution.getBinSize())
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
            unsigned long currentContigSize = 0;
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
                currentContigSize += count;
                back().push_back(currentBinIndex);
            }
            // bins do not spread across contigs
            if (!skipEmptyBins || currentContigSize)
            {
                // we don't create bins for empty contigs.
                ++currentBinIndex;
            }
        }
    }

    /**
     ** \brief convert a reference position on a contig into a bin index that
     ** can be used to identify either the file path or the stream associated
     ** to the ReferencePosition.
     **/
    size_t getBinIndex(const isaac::reference::ReferencePosition &referencePosition) const
    {
        const unsigned long contigId = referencePosition.getContigId();
        const std::vector<unsigned> &binIndexList = at(contigId + 1);
        const unsigned long position = referencePosition.getPosition();
        const unsigned long index = position / distributionBinSize_;
        assert(binIndexList.size() > index);
        return binIndexList[index];
    }

    /**
     * \return The first reference position that can be found in the bin
     */
    isaac::reference::ReferencePosition getBinFirstPos(const unsigned bin) const
    {
        std::vector<unsigned>::const_reference (std::vector<unsigned>::*f)() const = &std::vector<unsigned>::front;

        // user upper_bound to skip all the contigs that were so empty that they did not get mapped to a bin
        const_iterator binContigIterator = std::upper_bound(begin(), end(), bin, boost::bind(f, _2) > _1);
        ISAAC_ASSERT_MSG(begin() != binContigIterator, "Bin number has to be one of those we have a contig for: " << bin);
        //take a step back as we just skipped the last one we were looking for
        --binContigIterator;

        std::vector<unsigned>::const_iterator contigBinIterator =
            lower_bound(binContigIterator->begin(), binContigIterator->end(), bin);

        ISAAC_ASSERT_MSG(binContigIterator->end() != contigBinIterator, "Bin number must be present in the contig bins");

        const unsigned long contigId = binContigIterator - begin() - 1;
        const unsigned long position = distributionBinSize_ * (contigBinIterator - binContigIterator->begin());

        return isaac::reference::ReferencePosition(contigId, position);
    }

    /**
     * \return The first reference position that belongs to the subsequent bin. NOTE: for last bin in the contig
     *         there is not guarantee that no alignments will exist at this position and beyond. However, the amount
     *         of data aligning there should be considered minor and belonging to the last bin.
     */
    isaac::reference::ReferencePosition getBinFirstInvalidPos(const unsigned bin) const
    {
        std::vector<unsigned>::const_reference (std::vector<unsigned>::*f)() const = &std::vector<unsigned>::front;

        // user upper_bound to skip all the contigs that were so empty that they did not get mapped to a bin
        const_iterator binContigIterator = std::upper_bound(begin(), end(), bin, boost::bind(f, _2) > _1);
        ISAAC_ASSERT_MSG(begin() != binContigIterator, "Bin number has to be one of those we have a contig for: " << bin);
        //take a step back as we just skipped the last one we were looking for
        --binContigIterator;

        std::vector<unsigned>::const_iterator contigNextBinIterator =
            upper_bound(binContigIterator->begin(), binContigIterator->end(), bin);

        const unsigned long contigId = binContigIterator - begin() - 1;
        const unsigned long position = distributionBinSize_ * (contigNextBinIterator - binContigIterator->begin());

        return isaac::reference::ReferencePosition(contigId, position);
    }

    /**
     * \return The highest bin index to which the mapping is stored. Notice that there might be no bin with
     *         this index as it could have had no matches.
     */
    unsigned getHighestBinIndex() const
    {
        return back().back();
    }

};


inline std::ostream& operator << (std::ostream& os, const BinIndexMap &binIndexMap)
{
    std::string message = "\n";
    unsigned index = 0;
    BOOST_FOREACH(const std::vector<unsigned> &contigIndexList, binIndexMap)
    {
        message += boost::lexical_cast<std::string>(index) + ":";
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
        ++index;
    }
    os << message;
    return os;
}

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BIN_INDEX_MAP_HH
