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
 ** Component containing the data associated to a cluster: sequence and quality
 ** strings for all the reads in the cluster.
 **
 ** \author Come Raczy
 **/

#include <boost/foreach.hpp>

#include "alignment/Cluster.hh"
#include "flowcell/ReadMetadata.hh"
#include "oligo/Nucleotides.hh"

namespace isaac
{
namespace alignment
{

Cluster::Cluster(const unsigned maxReadLength)
    : tile_(0)
    , id_(0)
    , pf_(false)
    , barcodeLength_(0)
    , nonEmptyReads_(0)
{
    push_back(Read(maxReadLength, 0));
    push_back(Read(maxReadLength, 1));
}

void Cluster::init(
    const flowcell::ReadMetadataList &readMetadataList,
    std::vector<char>::const_iterator bclData,
    const unsigned tile,
    const unsigned long id,
    const ClusterXy &xy,
    const bool pf,
    const int barcodeLength)
{
    tile_ = tile;
    id_ = id;
    xy_ = xy;
    pf_ = pf;
    barcodeLength_ = barcodeLength;
    nonEmptyReads_ = 0;
    bclData_ = bclData;
    BOOST_FOREACH(const flowcell::ReadMetadata &readMetadata, readMetadataList)
    {
        //TODO: this check is to potentially allow 0-lenght reads (such as masked out) in readMetadataList
        // however, at the moment masked-out reads don't make into readMetadatList.
        if (readMetadata.getLength())
        {
            at(readMetadata.getIndex()).decodeBcl(bclData + barcodeLength_, bclData + readMetadata.getLength() + barcodeLength_, readMetadata.getIndex());
            bclData += readMetadata.getLength();
            ++nonEmptyReads_;
        }
    }
}

std::vector<char>::const_iterator Cluster::getBclData(const unsigned readIndex) const
{
    std::size_t offset = 0;
    for (unsigned i = 0; readIndex > i; ++i)
    {
        offset += (*this)[i].getForwardSequence().size();
    }
    return bclData_ + barcodeLength_ + offset;
}

unsigned long Cluster::getBarcodeSequence() const
{
    return oligo::packBclBases(bclData_, bclData_ + barcodeLength_);
}

} // namespace alignment
} // namespace isaac
