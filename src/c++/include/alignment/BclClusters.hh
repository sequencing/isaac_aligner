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
 ** \file BclClusters.hh
 **
 ** \brief In-memory representation of sequencing data
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_BCL_CLUSTERS_HH
#define iSAAC_ALIGNMENT_BCL_CLUSTERS_HH

#include <vector>


namespace isaac
{
namespace alignment
{

class BclClusters : std::vector<char>
{
    std::size_t clusterLength_;
    std::vector<bool> pf_;
public:
    using std::vector<char>::iterator;
    using std::vector<char>::const_iterator;
    BclClusters(unsigned clusterLength):
        clusterLength_(clusterLength)
    {
    }

    void reserveClusters(std::size_t reserveClusters)
    {
        reserve(clusterLength_ * reserveClusters);
        pf_.reserve(reserveClusters);
    }

    unsigned getClusterCount() const
    {
        return size() / clusterLength_;
    }

    unsigned getClusterLength() const
    {
        return clusterLength_;
    }

    /**
     * \post If the size of the buffer reduces, the data already in the buffer stays there.
     **/
    void reset(const std::size_t clusterLength, const std::size_t clusters)
    {
        clusterLength_ = clusterLength;
        resize(clusterLength_ * clusters);
        pf_.resize(clusters);
    }

    using std::vector<char>::end;
    iterator cluster(std::size_t cluster)
    {
        return begin() + cluster * clusterLength_;
    }
    const_iterator cluster(std::size_t cluster) const
    {
        return begin() + cluster * clusterLength_;
    }

    std::vector<bool> &pf()
    {
        return pf_;
    }

    bool pf(std::size_t cluster) const
    {
        return pf_.at(cluster);
    }
};

} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_BCL_CLUSTERS_HH
