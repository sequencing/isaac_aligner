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
 ** \file TileBarcodeStats.hh
 **
 ** \brief Statistics collection helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_ALIGNMENT_MATCH_SELECTOR_TILE_BARCODE_STATS_H
#define ISAAC_ALIGNMENT_MATCH_SELECTOR_TILE_BARCODE_STATS_H

#include "alignment/TemplateLengthStatistics.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

enum TemplateAlignmentType
{
    Normal,
    NmNm,       // seeds have no match in the reference.
    Qc,         // all seeds contain Ns, alignment is not possible
    Rm,         // one of the seeds exactly mapped to a high repeat or too many neighbors with the same prefix
    Filtered,   // user parameters required the template to be excluded
};

struct TileBarcodeStats
{
    explicit TileBarcodeStats() : templateLengthStatistics_()
    {
        reset();
    }

    void reset()
    {
        yield_ = 0;
        yieldQ30_ = 0;
        qualityScoreSum_ = 0;
        clusterCount_ = 0;
        unanchoredClusterCount_ = 0;
        nmnmClusterCount_ = 0;
        rmClusterCount_ = 0;
        qcClusterCount_ = 0;
        uniquelyAlignedFragmentCount_ = 0;
        alignedFragmentCount_ = 0;
        uniquelyAlignedPerfectFragmentCount_ = 0;
        alignmentScoreSum_ = 0;
        basesOutsideIndels_ = 0;
        uniquelyAlignedBasesOutsideIndels_ = 0;
        mismatches_ = 0;
        uniquelyAlignedMismatches_ = 0;
        fragmentCount_ = 0;
        templateLengthStatistics_.clear();
        templateLengthStatisticsSet_ = false;
        templateLengthStatisticsConflicts_ = false;
        alignmentModelCounts_[TemplateLengthStatistics::FFp] = 0;
        alignmentModelCounts_[TemplateLengthStatistics::FRp] = 0;
        alignmentModelCounts_[TemplateLengthStatistics::RFp] = 0;
        alignmentModelCounts_[TemplateLengthStatistics::RRp] = 0;
        alignmentModelCounts_[TemplateLengthStatistics::FFm] = 0;
        alignmentModelCounts_[TemplateLengthStatistics::FRm] = 0;
        alignmentModelCounts_[TemplateLengthStatistics::RFm] = 0;
        alignmentModelCounts_[TemplateLengthStatistics::RRm] = 0;

        std::fill(nominalModelCounts_, nominalModelCounts_ + TemplateLengthStatistics::CheckModelLast, 0);
    }

    unsigned long yield_;
    unsigned long yieldQ30_;
    unsigned long qualityScoreSum_;
    unsigned long clusterCount_;
    unsigned long unanchoredClusterCount_;
    unsigned long nmnmClusterCount_;
    unsigned long rmClusterCount_;
    unsigned long qcClusterCount_;
    unsigned long alignedFragmentCount_;
    unsigned long uniquelyAlignedFragmentCount_;
    // number of uniquely aligned fragments that have edit distance of 0
    unsigned long uniquelyAlignedPerfectFragmentCount_;
    unsigned long alignmentScoreSum_;
    unsigned long basesOutsideIndels_;
    unsigned long uniquelyAlignedBasesOutsideIndels_;
    unsigned long mismatches_;
    unsigned long uniquelyAlignedMismatches_;
    unsigned long alignmentModelCounts_[TemplateLengthStatistics::InvalidAlignmentModel+1];
    unsigned long nominalModelCounts_[TemplateLengthStatistics::CheckModelLast];

    unsigned long fragmentCount_;

    TemplateLengthStatistics templateLengthStatistics_;
    /**
     * indicates that this structure had the template length statistics set
     * (either via explicit recordTemplateLengthStatistics call or during the aggregation)
     */
    bool templateLengthStatisticsSet_;
    /**
     * indicates that the structure is an aggregation of multiple barcodes that have a conflicting template lengh stats
     */
    bool templateLengthStatisticsConflicts_;

    template<typename TemplateT>
    void recordTemplate(const TemplateT &templ)
    {
        ++alignmentModelCounts_[templ.getAlignmentModel()];
        const TemplateLengthStatistics::CheckModelResult cm = templ.checkModel();
        ++nominalModelCounts_[cm];
        ++clusterCount_;
        unanchoredClusterCount_ += templ.isUnanchored();
        nmnmClusterCount_ += templ.isNmNm();
        rmClusterCount_ += templ.isRm();
        qcClusterCount_ += templ.isQc();
    }
    template<typename FragmentT>
    void recordFragment(const FragmentT &fragment, const flowcell::ReadMetadata &readMetadata)
    {
        ISAAC_ASSERT_MSG(readMetadata.getLength() == fragment.getYield(), "Expected read length to match the read yield");
        yield_ += fragment.getYield();
        yieldQ30_ += fragment.getYieldQ30();
        qualityScoreSum_ += fragment.getQualityScoreSum();
        ++fragmentCount_;
        if (fragment.isAligned())
        {
            if (fragment.hasAlignmentScore())
            {
                alignmentScoreSum_ += fragment.getAlignmentScore();
            }
            mismatches_ += fragment.getMismatches();
            basesOutsideIndels_ += fragment.getUniquelyAlignedBasesOutsideIndels();
            ++alignedFragmentCount_;
        }

        if (fragment.isUniquelyAligned())
        {
            uniquelyAlignedMismatches_ += fragment.getMismatches();
            ++uniquelyAlignedFragmentCount_;
            uniquelyAlignedBasesOutsideIndels_ += fragment.getUniquelyAlignedBasesOutsideIndels();
            if (!fragment.getEditDistance())
            {
                ++uniquelyAlignedPerfectFragmentCount_;
            }
        }

    }

    void recordTemplateLengthStatistics(const TemplateLengthStatistics &templateLengthStatistics)
    {
        ISAAC_ASSERT_MSG(!templateLengthStatisticsSet_, "Setting template length stats is expected to happen only once");
        templateLengthStatistics_ = TemplateLengthStatistics(
            templateLengthStatistics.getMin(),
            templateLengthStatistics.getMax(),
            templateLengthStatistics.getMedian(),
            templateLengthStatistics.getLowStdDev(),
            templateLengthStatistics.getHighStdDev(),
            templateLengthStatistics.getBestModel(0),
            templateLengthStatistics.getBestModel(1),
            templateLengthStatistics.isStable());
        templateLengthStatisticsSet_ = true;
    }

    const TileBarcodeStats &operator +=(const TileBarcodeStats &right)
    {
        alignmentModelCounts_[TemplateLengthStatistics::FFp] += right.alignmentModelCounts_[TemplateLengthStatistics::FFp];
        alignmentModelCounts_[TemplateLengthStatistics::FRp] += right.alignmentModelCounts_[TemplateLengthStatistics::FRp];
        alignmentModelCounts_[TemplateLengthStatistics::RFp] += right.alignmentModelCounts_[TemplateLengthStatistics::RFp];
        alignmentModelCounts_[TemplateLengthStatistics::RRp] += right.alignmentModelCounts_[TemplateLengthStatistics::RRp];
        alignmentModelCounts_[TemplateLengthStatistics::FFm] += right.alignmentModelCounts_[TemplateLengthStatistics::FFm];
        alignmentModelCounts_[TemplateLengthStatistics::FRm] += right.alignmentModelCounts_[TemplateLengthStatistics::FRm];
        alignmentModelCounts_[TemplateLengthStatistics::RFm] += right.alignmentModelCounts_[TemplateLengthStatistics::RFm];
        alignmentModelCounts_[TemplateLengthStatistics::RRm] += right.alignmentModelCounts_[TemplateLengthStatistics::RRm];

        nominalModelCounts_[TemplateLengthStatistics::Oversized] += right.nominalModelCounts_[TemplateLengthStatistics::Oversized];
        nominalModelCounts_[TemplateLengthStatistics::Undersized] += right.nominalModelCounts_[TemplateLengthStatistics::Undersized];
        nominalModelCounts_[TemplateLengthStatistics::Nominal] += right.nominalModelCounts_[TemplateLengthStatistics::Nominal];
        nominalModelCounts_[TemplateLengthStatistics::NoMatch] += right.nominalModelCounts_[TemplateLengthStatistics::NoMatch];

        yield_ += right.yield_;
        yieldQ30_ += right.yieldQ30_;
        qualityScoreSum_ += right.qualityScoreSum_;
        clusterCount_ += right.clusterCount_;
        unanchoredClusterCount_ += right.unanchoredClusterCount_;
        nmnmClusterCount_ += right.nmnmClusterCount_;
        rmClusterCount_ += right.rmClusterCount_;
        qcClusterCount_ += right.qcClusterCount_;
        uniquelyAlignedFragmentCount_ += right.uniquelyAlignedFragmentCount_;
        alignedFragmentCount_ += right.alignedFragmentCount_;
        uniquelyAlignedPerfectFragmentCount_ += right.uniquelyAlignedPerfectFragmentCount_;
        alignmentScoreSum_ += right.alignmentScoreSum_;
        basesOutsideIndels_ += right.basesOutsideIndels_;
        uniquelyAlignedBasesOutsideIndels_ += right.uniquelyAlignedBasesOutsideIndels_;
        mismatches_ += right.mismatches_;
        uniquelyAlignedMismatches_ += right.uniquelyAlignedMismatches_;
        fragmentCount_ += right.fragmentCount_;

        if (!templateLengthStatisticsConflicts_)
        {
            // these are not accumulated. They are propagated.
            if (!templateLengthStatisticsSet_)
            {
                templateLengthStatistics_ = right.templateLengthStatistics_;
                templateLengthStatisticsSet_ = right.templateLengthStatisticsSet_;
            }
            else
            {
                //At the moment, there is no valid scenario in which two stable statistics get aggregated
                templateLengthStatisticsConflicts_ = right.templateLengthStatisticsSet_;
            }
        }
        return *this;
    }

    void finalize()
    {
    }

};

} //namespace matchSelector
} //namespace alignment
} //namespace isaac

#endif //ISAAC_ALIGNMENT_MATCH_SELECTOR_TILE_BARCODE_STATS_H
