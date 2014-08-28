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
 ** \file MatchDistribution.hh
 **
 ** \brief Tracking of the distribution of the matches across the reference
 ** genome.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_DISTRIBUTION_HH
#define iSAAC_ALIGNMENT_MATCH_DISTRIBUTION_HH

#include <vector>

#include "common/Debug.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Counts of matches for bins of regular size along the genome.
 **
 ** The external vector stores one vector for each contig in the reference
 ** genome. Each internal 'std::vector<unsigned>' represents the count of
 ** matches for all the bins in the corresponding contig.
 **
 ** The bin sizes is as determined by MatchDistribution::binSize.
 **/
class MatchDistribution: public std::vector<std::vector<unsigned> >
{
public:
    MatchDistribution()
    {

    }
    MatchDistribution(const isaac::reference::SortedReferenceMetadataList &sortedReferenceMetadataList)
    {
        initialize(sortedReferenceMetadataList);
    }

    /**
     ** \brief Initialize the MatchDistribution for the specified references.
     **
     ** Creates one vector for each contig, using the max contig length of each contig
     ** from all the provided references to infer the number of bins,
     ** and initialize all bin counts to 0.
     **/
    void initialize(const isaac::reference::SortedReferenceMetadataList &sortedReferenceMetadataList);

    /**
     ** \brief Consolidate all the partial match distribution
     ** given as a parameter.
     **
     ** simply add the values from the partial distributions of corresponding contigs.
     **/
    void consolidate(const MatchDistribution &matchDistribution);

    /**
     ** \brief Consolidate all the partial match distributions from the list
     ** given as a parameter.
     **/
    void consolidate(const std::vector<MatchDistribution> &matchDistributionList);

    /**
     ** \brief Add a match for the indicated (contig, position)
     **/
    void addMatches(const size_t contigIndex, const size_t position, const unsigned count)
    {
        ISAAC_ASSERT_MSG(size() > contigIndex, "Contig index too large");
        const size_t binIndex = getBinIndex(position);
        ISAAC_ASSERT_MSG((*this)[contigIndex].size() > binIndex, "Contig bin has not been initialized to handle the position");
        (*this)[contigIndex][binIndex] += count;
    }

    /**
     * Bins must be granular enough to allow MatchSelector binning flexibility for highly-covered tiny genomes (e.g. PhiX).
     * This presents a bit of a memory challenge when finding matches for large genomes on big number of threads.
     * For example:
     *  64 cores, 2^32 bases long human genome, 2^8 Bin size:
     *      64 * 2^32 / 2^8 * sizeof(unsigned) = 4 gigabytes. And we need those gigabytes for seeds, reference, etc...
     * So, keep it 2^11 for now
     */
    unsigned long getBinSize() const {return (1 << 11);}

    bool isEmptyContig(const size_t contigIndex) const
    {
        const std::vector<unsigned> &contigBins = (*this)[contigIndex];
        return contigBins.end() == std::find_if(contigBins.begin(), contigBins.end(), boost::bind(std::not_equal_to<unsigned>(), _1, 0));
    }
private:

    /// convert a position on a contig into a bin index
    size_t getBinIndex(unsigned long position) const {return position / getBinSize();}
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_DISTRIBUTION_HH
