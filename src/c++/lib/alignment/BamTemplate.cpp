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
 ** \file BamTemplate.cpp
 **
 ** \brief See BamTemplate.hh
 ** 
 ** \author Come Raczy
 **/

#include <boost/foreach.hpp>

#include "alignment/BamTemplate.hh"

namespace isaac
{
namespace alignment
{

BamTemplate::BamTemplate(const std::vector<unsigned> &cigarBuffer)
    : cigarBuffer_(cigarBuffer)
    , alignmentScore_(0)
{
    fragmentMetadataList_.reserve(2);
}

void BamTemplate::initialize(const flowcell::ReadMetadataList &tileReads, const Cluster &cluster)
{
    fragmentMetadataList_.clear();
    alignmentScore_ = 0;
    properPair_ = false;
    BOOST_FOREACH(const flowcell::ReadMetadata &read, tileReads)
    {
        fragmentMetadataList_.push_back(FragmentMetadata(&cluster, &cigarBuffer_, read.getIndex()));
    }
}

bool BamTemplate::filterLowQualityFragments(const unsigned mapqThreshold)
{
    bool ret = false;
    unsigned alignmentScore = 0;
    for (unsigned i = 0; getFragmentCount() > i; ++i)
    {
        FragmentMetadata &fragment = getFragmentMetadata(i);
        if (mapqThreshold > fragment.getAlignmentScore())
        {
            fragment.cigarLength = 0;
            fragment.cigarOffset = 0;
            fragment.alignmentScore = 0;
            const FragmentMetadata &mate = getFragmentMetadata((i+1) % getFragmentCount());
            fragment.position = mate.position;
            fragment.contigId = mate.contigId;
        }
        else if (fragment.isAligned())
        {
            ret = true;
        }
        // update the alignment score of the templates that didn't resolve into a proper pair
        alignmentScore += fragment.alignmentScore;
    }
    setAlignmentScore(alignmentScore);
    return ret;
}

} // namespace alignment
} // namespace isaac
