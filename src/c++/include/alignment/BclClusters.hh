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
 ** \file BclClusters.hh
 **
 ** \brief In-memory representation of sequencing data
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_BCL_CLUSTERS_HH
#define iSAAC_ALIGNMENT_BCL_CLUSTERS_HH

#include <vector>

#include "alignment/Cluster.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

class BclClusters : std::vector<char>
{
    std::size_t clusterLength_;
    std::vector<bool> pf_;
    std::vector<ClusterXy> xy_;
public:
    using std::vector<char>::iterator;
    using std::vector<char>::const_iterator;
    BclClusters(unsigned clusterLength):
        clusterLength_(clusterLength)
    {
    }

    void reserveClusters(const std::size_t reserveClusters, const bool storeXy)
    {
        reserve(clusterLength_ * reserveClusters);
        pf_.reserve(reserveClusters);
        if (storeXy)
        {
            xy_.reserve(reserveClusters);
        }
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
        if (storeXy())
        {
            xy_.resize(clusters);
        }
    }

    void reduceWastedMemory()
    {
        ISAAC_THREAD_CERR << "BclClusters reducing memory waste" << std::endl;
        std::vector<char> bclTmp(*this);
        swap(bclTmp);
        std::vector<bool> pfTmp(pf_);
        pf_.swap(pfTmp);
        ISAAC_THREAD_CERR << "BclClusters reducing memory waste done. Saved: " <<
            bclTmp.capacity() - capacity() + pfTmp.capacity() / 8 - pf_.capacity() / 8 << " bytes" << std::endl;

        // no reducing memory on xy_ as it is not supposed to be used in this situation
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

    std::vector<ClusterXy> &xy()
    {
        return xy_;
    }

    const ClusterXy &xy(const std::size_t cluster) const
    {
        static const ClusterXy unsetXy;
        return storeXy() ? xy_.at(cluster) : unsetXy;
    }

    bool storeXy() const {return xy_.capacity();}
};

} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_BCL_CLUSTERS_HH
