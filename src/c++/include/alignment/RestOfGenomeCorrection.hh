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
 ** \file RestOfGenomeCorrection.hh
 **
 ** \brief Component to encapsulate probability of read randomly aligning to the genome.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_REST_OF_GENOME_CORRECTION_HH
#define iSAAC_ALIGNMENT_REST_OF_GENOME_CORRECTION_HH

#include <boost/foreach.hpp>

#include "common/Debug.hh"
#include "reference/Contig.hh"

namespace isaac
{
namespace alignment
{


class RestOfGenomeCorrection
{
public:
    RestOfGenomeCorrection(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList)
    {
        setGenome(contigList, readMetadataList);
    }

    void setGenome(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList)
    {
        getRestOfGenomeCorrection(contigList, readMetadataList, rogCorrectionList_);
        rogCorrection_ = Quality::restOfGenomeCorrection(reference::genomeLength(contigList),
                                                         flowcell::getTotalReadLength(readMetadataList));
        // can't have 0.0 in it as it turns the alignment score into 0
        rogCorrection_ = std::max(rogCorrection_, std::numeric_limits<double>::min());
    }

    double getReadRogCorrection(const unsigned readIndex) const
    {
        ISAAC_ASSERT_MSG(sizeof(rogCorrectionList_) > readIndex, "Only up to 2 reads supported");
        return rogCorrectionList_[readIndex];
    }
    double getRogCorrection() const {return rogCorrection_;}
private:
    static const unsigned READS_MAX = 2;
    /// Rest-of-genome correction for individual fragments
    double rogCorrectionList_[READS_MAX];
    /// Rest-of-genome correction for the template when all fragments match
    double rogCorrection_;


    static void getRestOfGenomeCorrection(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        double rogCorrectionList[READS_MAX])
    {
        const size_t genomeLength = reference::genomeLength(contigList);
        BOOST_FOREACH(const flowcell::ReadMetadata &readMetadata, readMetadataList)
        {
            rogCorrectionList[readMetadata.getIndex()] = Quality::restOfGenomeCorrection(genomeLength, readMetadata.getLength());
            // can't have 0.0 in it as it turns the alignment score into 0
            rogCorrectionList[readMetadata.getIndex()] = std::max(rogCorrectionList[readMetadata.getIndex()], std::numeric_limits<double>::min());
        }
    }
};


} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_REST_OF_GENOME_CORRECTION_HH
