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
 ** \file Cluster.hh
 **
 ** \brief Component containing the data associated to a cluster: sequence and
 ** quality strings for all the reads in the cluster.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_CLUSTER_HH
#define iSAAC_ALIGNMENT_CLUSTER_HH

#include <vector>
#include <string>

#include "alignment/Read.hh"
#include "common/FastIo.hh"

namespace isaac
{

namespace flowcell
{

class ReadMetadataList;

}

namespace alignment
{
struct ClusterXy
{
    ClusterXy() : x_(POSITION_NOT_SET), y_(POSITION_NOT_SET){}
    ClusterXy(int x, int y) : x_(x), y_(y) {}
    ClusterXy(const std::pair<int, int> &xy) : x_(xy.first), y_(xy.second) {}
    static const int POSITION_NOT_SET = boost::integer_traits<int>::const_max;

    int x_;
    int y_;

    bool isSet() const {return POSITION_NOT_SET != x_;}
};

class Cluster: public std::vector<Read>
{
public:
    Cluster(const unsigned maxReadLen);

    void init(
        const flowcell::ReadMetadataList &readMetadataList,
        std::vector<char>::const_iterator bclData,
        const unsigned tile,
        const unsigned long id,
        const ClusterXy &xy,
        const bool pf,
        const int barcodeLength);
    unsigned long getTile() const {return tile_;}
    unsigned long getId() const {return id_;}
    bool getPf() const {return pf_;}
    ClusterXy getXy() const {return xy_;}
    int getBarcodeLength() const { return barcodeLength_; }
    unsigned getNonEmptyReadsCount() const {return nonEmptyReads_;}
    /// Beginning of the BCL data for the indicated read
    std::vector<char>::const_iterator getBclData(unsigned readIndex) const;
    unsigned long getBarcodeSequence() const;
    template<class InpuT> friend InpuT& operator >>(InpuT &input, Cluster &cluster);
private:
    unsigned tile_;
    unsigned long id_;
    ClusterXy xy_;
    bool pf_;
    int barcodeLength_;
    unsigned nonEmptyReads_;
    std::vector<char>::const_iterator bclData_;
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_CLUSTER_HH
