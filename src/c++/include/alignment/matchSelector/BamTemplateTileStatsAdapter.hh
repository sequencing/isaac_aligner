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
 ** \file BamTemplateTileStatsAdapter.hh
 **
 ** \brief Conversion from BamTemplate to the interface suitable for MatchSelectorTileStats::recordTemplate
 ** generation.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BAM_TEMPLATE_TILE_STATS_ADAPTER_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_BAM_TEMPLATE_TILE_STATS_ADAPTER_HH

#include <boost/noncopyable.hpp>

#include "alignment/BamTemplate.hh"
#include "alignment/matchSelector/TileBarcodeStats.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

class BamTemplateTileStatsAdapter: boost::noncopyable
{
    const TemplateLengthStatistics &templateLengthStatistics_;
    const BamTemplate &templ_;
    const TemplateAlignmentType templateType_;
public:
    BamTemplateTileStatsAdapter(const TemplateLengthStatistics &templateLengthStatistics,
                                const BamTemplate &templ,
                                const TemplateAlignmentType templateType):
        templateLengthStatistics_(templateLengthStatistics), templ_(templ), templateType_(templateType){}

    TemplateLengthStatistics::AlignmentModel getAlignmentModel() const
    {
        ISAAC_ASSERT_MSG(2 >= templ_.getFragmentCount(), "single ended or paired data supported only");
        if (2 == templ_.getFragmentCount())
        {
            // compatibility with kagu statistics: only pairs of uniquely aligned fragments count
            return (templ_.getFragmentMetadata(0).isUniquelyAligned() &&
                    templ_.getFragmentMetadata(1).isUniquelyAligned())
                    ? TemplateLengthStatistics::alignmentModel(templ_.getFragmentMetadata(0), templ_.getFragmentMetadata(1))
                    : TemplateLengthStatistics::InvalidAlignmentModel;
        }
        else
        {
            return TemplateLengthStatistics::InvalidAlignmentModel;
        }
    }

    TemplateLengthStatistics::CheckModelResult checkModel() const
    {
        ISAAC_ASSERT_MSG(2 >= templ_.getFragmentCount(), "single ended or paired data supported only");
        if (2 == templ_.getFragmentCount())
        {
            // compatibility with kagu statistics: only pairs of uniquely aligned fragments count
            return (templ_.getFragmentMetadata(0).isUniquelyAligned() &&
                    templ_.getFragmentMetadata(1).isUniquelyAligned())
                    ? templateLengthStatistics_.checkModel(templ_.getFragmentMetadata(0), templ_.getFragmentMetadata(1)):
                    TemplateLengthStatistics::NoMatch;
        }
        else
        {
            return TemplateLengthStatistics::NoMatch;
        }
    }

    unsigned long getMismatches() const
    {
        return templ_.getMismatchCount();
    }

    unsigned getAlignmentScore() const
    {
        return templ_.getAlignmentScore();
    }

    bool hasAlignmentScore() const
    {
        return templ_.hasAlignmentScore();
    }

    bool isUnanchored() const
    {
        return templ_.isUnanchored();
    }

    bool isNmNm() const
    {
        return NmNm == templateType_;
    }

    bool isRm() const
    {
        return  Rm == templateType_;
    }

    bool isQc() const
    {
        return  Qc == templateType_;
    }
};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BAM_TEMPLATE_TILE_STATS_ADAPTER_HH
