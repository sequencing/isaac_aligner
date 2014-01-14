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
 ** \file MatchSelectorStats.hh
 **
 ** \brief FragmentDispatcher statistics helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_DISPATCHER_STATS_H
#define ISAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_DISPATCHER_STATS_H

#include "alignment/BamTemplate.hh"
#include "alignment/matchSelector/TileStats.hh"
#include "alignment/matchSelector/TileBarcodeStats.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "flowcell/TileMetadata.hh"

#include "FragmentMetadataTileStatsAdapter.hh"
#include "BamTemplateTileStatsAdapter.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

namespace fs=boost::filesystem;

class MatchSelectorStats
{
public:

    MatchSelectorStats(
        const flowcell::BarcodeMetadataList &barcodeMetadataList) :
            barcodeMetadataList_(barcodeMetadataList)
    {
        ISAAC_TRACE_STAT("MatchSelectorStats::MatchSelectorStats ")
        const unsigned tileStatsCount = maxReads_ * filterStates_;
        ISAAC_THREAD_CERR << "Allocating " << tileStatsCount << " tile stats." << std::endl;
        tileStats_.resize(tileStatsCount);
        ISAAC_THREAD_CERR << "Allocating " << tileStatsCount << " tile stats done. Total size is " << tileStats_.capacity() * sizeof(TileStats) << " bytes."<< std::endl;

        // TODO: this allows for a situation where every barcode is expected to be found on every tile.
        // Since barcodes are constrained to one lane, plenty of ram can be saved by coming up with an alternative solution.
        const unsigned tileBarcodeStatsCount = maxReads_ * barcodeMetadataList.size() * filterStates_;
        ISAAC_THREAD_CERR << "Allocating " << tileBarcodeStatsCount << " tile barcode stats." << std::endl;
        tileBarcodeStats_.resize(tileBarcodeStatsCount);
        ISAAC_THREAD_CERR << "Allocating " << tileBarcodeStatsCount << " tile barcode stats done. Total size is " << tileBarcodeStats_.capacity() * sizeof(TileBarcodeStats) << " bytes."<< std::endl;
        ISAAC_TRACE_STAT("MatchSelectorStats::MatchSelectorStats constructed")
    }

    void reset()
    {
        std::for_each(tileStats_.begin(), tileStats_.end(),
                      boost::bind(&TileStats::reset, _1));
        std::for_each(tileBarcodeStats_.begin(), tileBarcodeStats_.end(),
                      boost::bind(&TileBarcodeStats::reset, _1));

    }

    void recordTemplate(
        const flowcell::ReadMetadataList &readMetadatalist,
        const TemplateLengthStatistics &templateLengthStatistics,
        const BamTemplate &bamTemplate,
        const unsigned barcodeIndex,
        const TemplateAlignmentType templateType)
    {
        // all pair-level stats are recorded under read index 0
        BamTemplateTileStatsAdapter tileStatsAdapter(templateLengthStatistics, bamTemplate, templateType);
        if (bamTemplate.getPassesFilter())
        {
            tileStats_.at(tileIndex(bamTemplate.getFragmentMetadata(0), true)).recordTemplate(tileStatsAdapter);
            tileBarcodeStats_.at(tileBarcodeIndex(bamTemplate.getFragmentMetadata(0), barcodeIndex, true)).recordTemplate(tileStatsAdapter);
        }
        tileStats_.at(tileIndex(bamTemplate.getFragmentMetadata(0), false)).recordTemplate(tileStatsAdapter);
        tileBarcodeStats_.at(tileBarcodeIndex(bamTemplate.getFragmentMetadata(0), barcodeIndex, false)).recordTemplate(tileStatsAdapter);
        for(unsigned i = 0; bamTemplate.getFragmentCount() > i; ++i)
        {
            const FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(i);
            FragmentMetadataTileStatsAdapter tileStatsAdapter(fragment);
            if (bamTemplate.getPassesFilter())
            {
                tileStats_.at(tileIndex(fragment, true)).recordFragment(tileStatsAdapter, readMetadatalist.at(i));
                tileBarcodeStats_.at(tileBarcodeIndex(fragment, barcodeIndex, true)).recordFragment(tileStatsAdapter, readMetadatalist.at(i));
            }
            tileStats_.at(tileIndex(fragment, false)).recordFragment(tileStatsAdapter, readMetadatalist.at(i));
            tileBarcodeStats_.at(tileBarcodeIndex(fragment, barcodeIndex, false)).recordFragment(tileStatsAdapter, readMetadatalist.at(i));
        }
    }

    void recordTemplateLengthStatistics(
        const flowcell::BarcodeMetadata &barcodeMetadata,
        const TemplateLengthStatistics &templateLengthStatistics)
    {
        tileBarcodeStats_.at(tileBarcodeIndex(barcodeMetadata)).recordTemplateLengthStatistics(templateLengthStatistics);
    }

    MatchSelectorStats &operator +=(const MatchSelectorStats &right)
    {
        ISAAC_ASSERT_MSG(right.barcodeMetadataList_.size() == barcodeMetadataList_.size(), "dimensions must match");
        ISAAC_ASSERT_MSG(right.tileBarcodeStats_.size() == tileBarcodeStats_.size(), "size must match");
        ISAAC_ASSERT_MSG(right.tileStats_.size() == tileStats_.size(), "size must match");
        std::transform(tileStats_.begin(), tileStats_.end(),
                       right.tileStats_.begin(), tileStats_.begin(), std::plus<TileStats>());
        unsigned i = 0;
        BOOST_FOREACH(TileBarcodeStats &tileBarcodeStats, tileBarcodeStats_)
        {
            tileBarcodeStats += right.tileBarcodeStats_.at(i);
            ++i;
        }
        return *this;
    }

    const MatchSelectorStats operator +(const MatchSelectorStats &right) const
    {
        MatchSelectorStats ret(*this);
        ret += right;
        return ret;
    }

    MatchSelectorStats & operator =(const MatchSelectorStats &that) {
        ISAAC_ASSERT_MSG(that.barcodeMetadataList_.size() == barcodeMetadataList_.size(), "dimensions must match");
        ISAAC_ASSERT_MSG(that.tileStats_.size() == tileStats_.size(), "size must match");
        ISAAC_ASSERT_MSG(that.tileBarcodeStats_.size() == tileBarcodeStats_.size(), "size must match");
        tileStats_ = that.tileStats_;
        tileBarcodeStats_ = that.tileBarcodeStats_;
        return *this;
    }

    const TileBarcodeStats &getReadBarcodeTileStat(
        const flowcell::ReadMetadata& read,
        const flowcell::BarcodeMetadata& barcode,
        const bool passesFilter) const
    {
        return tileBarcodeStats_.at(tileBarcodeIndex(read, barcode, passesFilter));
    }

    const TileStats &getReadTileStat(
        const flowcell::ReadMetadata& read,
        const bool passesFilter) const
    {
        return tileStats_.at(tileIndex(read, passesFilter));
    }

    void finalize()
    {
        std::for_each(tileStats_.begin(), tileStats_.end(), boost::bind(&TileStats::finalize, _1));
        std::for_each(tileBarcodeStats_.begin(), tileBarcodeStats_.end(), boost::bind(&TileBarcodeStats::finalize, _1));
    }

private:
    static const unsigned filterStates_ = 2;
    static const unsigned maxReads_ = 2;
    const std::vector<flowcell::BarcodeMetadata> &barcodeMetadataList_;

    /**
     * \brief detailed blank, mismatch and alignment score distribution stats
     */
    std::vector<TileStats>  tileStats_;
    /**
     * \brief higher-level stats that we can afford to keep per tile-barcode
     */
    std::vector<TileBarcodeStats>  tileBarcodeStats_;

    unsigned tileBarcodeIndex(
        const flowcell::ReadMetadata& read,
        const flowcell::BarcodeMetadata& barcode,
        const bool passesFilter) const
    {
        return tileBarcodeIndex(read.getIndex(), barcode.getIndex(), passesFilter);
    }

    unsigned tileIndex(
        const flowcell::ReadMetadata& read,
        const bool passesFilter) const
    {
        return tileIndex(read.getIndex(), passesFilter);
    }

    unsigned tileBarcodeIndex(const flowcell::BarcodeMetadata& barcode) const
    {
        return tileBarcodeIndex(0, barcode.getIndex(), false);
    }

    unsigned tileBarcodeIndex(
        const FragmentMetadata &fragment,
        const unsigned barcodeIndex,
        bool passesFilter) const
    {
        return tileBarcodeIndex(fragment.getReadIndex(), barcodeIndex, passesFilter);
    }

    unsigned tileIndex(
        const FragmentMetadata &fragment,
        bool passesFilter) const
    {
        return tileIndex(fragment.getReadIndex(), passesFilter);
    }

    unsigned tileBarcodeIndex(const unsigned read, const unsigned barcode, const bool passesFilter) const
    {
        return
            read * barcodeMetadataList_.size() * filterStates_ +
            barcode * filterStates_ +
            passesFilter;
    }

    unsigned tileIndex(const unsigned read, const bool passesFilter) const
    {
        return
            read * filterStates_ +
            passesFilter;
    }
};

} //namespace matchSelector
} //namespace alignment
} //namespace isaac

#endif //ISAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_DISPATCHER_STATS_H
