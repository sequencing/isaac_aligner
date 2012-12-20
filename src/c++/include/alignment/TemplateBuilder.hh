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
 ** \file TemplateBuilder.hh
 **
 ** \brief Construction of BamTemplate instances
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_HH
#define iSAAC_ALIGNMENT_TEMPLATE_BUILDER_HH

#include <boost/noncopyable.hpp>

#include "reference/Contig.hh"
#include "flowcell/ReadMetadata.hh"
#include "alignment/SeedMetadata.hh"
#include "alignment/Match.hh"
#include "alignment/Cluster.hh"
#include "alignment/BamTemplate.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "alignment/BandedSmithWaterman.hh"
#include "alignment/ShadowAligner.hh"
#include "alignment/Cigar.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Utility component creating Template instances from Seed Matches.
 **
 ** The intended use is to create an instance of a TemplateBuilder for each
 ** thread and delegate to that instance the identification of the most likely
 ** template for each cluster. This is done by invoking the build method on the
 ** complete set of matches identified so far for the cluster. In the build
 ** method, the TemplateBuilder will do the alignment (first gapped, then
 ** ungapped), calculate the alignment quality of the individual fragments,
 ** select the most likely combination of fragments, resolve repeats and try
 ** aligning orphans.
 **
 ** TODO: add Template statistics
 **/
class TemplateBuilder: boost::noncopyable
{
public:
    enum DodgyAlignmentScore
    {
        Zero,
        Unknown,
        Unaligned
    };


    /**
     ** \brief Construct a template builder for a reference genome and a given
     ** set of reads.
     **
     ** TODO: Add support for clipping (variable length genomes)
     **/
    TemplateBuilder(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned repeatThreshold,
        const unsigned maxSeedsPerRead,
        const unsigned gappedMismatchesMax,
        const bool scatterRepeats,
        const int gapMatchScore = 2,
        const int gapMismatchScore = -1,
        const int gapOpenScore = -15,
        const int gapExtendScore = -3,
        const DodgyAlignmentScore dodgyAlignmentScore = Unaligned);

    bool buildFragments(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        const SeedMetadataList &seedMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<Match>::const_iterator matchBegin,
        const std::vector<Match>::const_iterator matchEnd,
        const Cluster &cluster,
        const bool withGaps)
    {
        return fragmentBuilder_.build(contigList, readMetadataList, seedMetadataList, sequencingAdapters,
                                       matchBegin, matchEnd, cluster, withGaps);
    }

    const std::vector<std::vector<FragmentMetadata> > &getFragments() const {return fragmentBuilder_.getFragments();}

    /**
     ** \brief Build the most likely template for a single cluster, givena set of fragments
     **
     ** \return false means template ended up not having a single read alignned anywhere.
     **
     ** This method will initialize the internal template of the builder.
     **
     ** Precondition: the input list of fragments is partitionned by readIndex
     ** and sorted genomic position. This means that the order is (tileId,
     ** clusterId, seedIndex, reverse, contig, position) where 'tileId' and
     ** 'clusterId' are constant.
     **/
    bool buildTemplate(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const Cluster &cluster,
        const TemplateLengthStatistics &templateLengthStatistics,
        const unsigned mapqThreshold);

    /**
     * \brief Same as above but unit testing friendly.
     */
    bool buildTemplate(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const Cluster &cluster,
        const TemplateLengthStatistics &templateLengthStatistics);

    /**
     ** \brief Getter for the BamTemplate
     **/
    const BamTemplate &getBamTemplate() const {return bamTemplate_;}
    BamTemplate &getBamTemplate() {return bamTemplate_;}
private:
    // when considering orphans for shadow alignment, don't look at those that are further than
    // orphanLogProbabilitySlack_ away from the best orphan
    static const double orphanLogProbabilitySlack_ = 100.0;
    static const unsigned readsMax_ = 2;
    // number of repeat alignment candidate locations that we are prepared to track. If it's above that
    // then the location will not be assigned
    static const unsigned TRACKED_REPEATS_MAX_EVER = 1000;
    static const unsigned SKIP_ORPHAN_EDIT_DISTANCE = 3;

    // Maximum alignment score given to fragments and templates that are not well anchored but don't have any mismatches
    static const unsigned DODGY_BUT_CLEAN_ALIGNMENT_SCORE = 10U;

    const bool scatterRepeats_;
    const DodgyAlignmentScore dodgyAlignmentScore_;

    /// Helper component to align fragments individually
    FragmentBuilder fragmentBuilder_;
    /// Cached storage for iterative template building
    BamTemplate bamTemplate_;
    /// Helper component to rescue shadow reads (or poorly aligned fragments)
    ShadowAligner shadowAligner_;
    /// Buffer for the cigar strings of rescued clusters
    std::vector<unsigned> cigarBuffer_;
    /// Buffer for the list of shadows rescued by the shadow aligner
    std::vector<FragmentMetadata> shadowList_;

    /// arrays temporary used in buildDisjointTemplate and rescueShadow

    struct PairProbability
    {
        reference::ReferencePosition r1Pos_;
        reference::ReferencePosition r2Pos_;
        double logProbability_;

        PairProbability() : r1Pos_(reference::ReferencePosition::NoMatch), r2Pos_(reference::ReferencePosition::NoMatch), logProbability_(0.0){}

        PairProbability(
                const reference::ReferencePosition r1Pos,
                const reference::ReferencePosition r2Pos,
                double logProbability) : r1Pos_(r1Pos), r2Pos_(r2Pos), logProbability_(logProbability){}

        bool operator < (const PairProbability &that) const
        {
            //if we happen to have pair and its inversion, let's have the higher probability on top so that it counts towards the total
            return r1Pos_ < that.r1Pos_ || (r1Pos_ == that.r1Pos_ && (r2Pos_ < that.r2Pos_ || (r2Pos_ == that.r2Pos_ && logProbability_ > that.logProbability_)));
        }
    };

    // Max number of orphans times the max number of shadows that can be rescued for each orphan times the number of reads
    common::FiniteCapacityVector<PairProbability, TRACKED_REPEATS_MAX_EVER * TRACKED_REPEATS_MAX_EVER * readsMax_> allPairProbabilities_;


    struct ShadowProbability
    {
        reference::ReferencePosition pos_;
        double logProbability_;

        ShadowProbability() : pos_(reference::ReferencePosition::NoMatch), logProbability_(0.0){}
        ShadowProbability(
                const reference::ReferencePosition pos,
                double logProbability) : pos_(pos), logProbability_(logProbability){}

        bool operator < (const ShadowProbability &that) const
        {
            //if we happen to have fragment and its inversion, let's have the higher probability on top so that it counts towards the total
            return pos_ < that.pos_ || (pos_ == that.pos_ && logProbability_ > that.logProbability_);
        }
    };

    // Max number of orphans times the max number of shadows that can be rescued for each orphan plus the max number of shadows that have been discovered during seed matching
    common::FiniteCapacityVector<ShadowProbability, TRACKED_REPEATS_MAX_EVER * TRACKED_REPEATS_MAX_EVER + TRACKED_REPEATS_MAX_EVER> allShadowProbabilities_[readsMax_];
    common::FiniteCapacityVector<FragmentMetadata, TRACKED_REPEATS_MAX_EVER> bestOrphanShadows_[readsMax_];

    /// Helper method to select the best fragment for single-ended runs
    bool pickBestFragment(
        const TemplateLengthStatistics &templateLengthStatistics,
        const std::vector<FragmentMetadata> &fragmentList);
    bool rescueShadow(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const TemplateLengthStatistics &templateLengthStatistics);
    /// Helper method to find the best pair of fragments for paired-end runs
    bool pickBestPair(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const TemplateLengthStatistics &templateLengthStatistics);

    typedef std::vector<FragmentMetadata>::const_iterator FragmentIterator;

    struct BestPairInfo
    {
        BestPairInfo(FragmentIterator bestR1Fragment, FragmentIterator bestR2Fragment) :
            bestTemplateLogProbability(-std::numeric_limits<double>::max()),
            bestTemplateScore(-1),
            resolvedTemplateCount(0), bestPairEditDistance(0),
            totalTemplateProbability(0.0)
        {
            bestPairFragments[0].push_back(bestR1Fragment);
            bestPairFragments[1].push_back(bestR2Fragment);
        }

        BestPairInfo() :
            bestTemplateLogProbability(-std::numeric_limits<double>::max()),
            bestTemplateScore(-1),
            resolvedTemplateCount(0), bestPairEditDistance(0),
            totalTemplateProbability(0.0)
        {
        }
        typedef isaac::common::FiniteCapacityVector<FragmentIterator, TRACKED_REPEATS_MAX_EVER> FragmentIteratorVector;
        FragmentIteratorVector bestPairFragments[readsMax_];
        double bestTemplateLogProbability;
        unsigned long bestTemplateScore;
        unsigned resolvedTemplateCount;
        unsigned bestPairEditDistance;
        double totalTemplateProbability;
    };
    friend std::ostream & operator << (std::ostream & os, const BestPairInfo& bestPairInfo);

    /// Helper method to locate the best pair of fragments
    BestPairInfo locateBestPair(
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const TemplateLengthStatistics &templateLengthStatistics) const;
    /// Helper method to build a paired-end template
    bool buildPairedEndTemplate(
        const TemplateLengthStatistics &templateLengthStatistics,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        BestPairInfo &bestPairInfo);
    /// Helper method to produce a template when no pair of fragments match
    bool buildDisjoinedTemplate(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const TemplateLengthStatistics &templateLengthStatistics,
        const BestPairInfo &knownBestPair);

    bool scoreDisjoinedTemplate(
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const TemplateLengthStatistics &templateLengthStatistics,
        const BestPairInfo &bestOrphans,
        const BestPairInfo &knownBestPair,
        const unsigned bestOrphanIndex,
        const double totalShadowProbability,
        const double totalOrphanProbability,
        const FragmentIterator bestDisjoinedFragments[readsMax_]);

    FragmentMetadata cloneWithCigar(const FragmentMetadata &right);

    /// Helper method to select the fragment with the highest logProbability
    std::vector<FragmentMetadata>::const_iterator getBestFragment(
        const std::vector<FragmentMetadata> &fragmentList) const;
    /// Helper function to calculate the alignment score of a fragment
    bool updateMappingScore(
        FragmentMetadata &fragment,
        const TemplateLengthStatistics &templateLengthStatistics,
        const FragmentIterator listFragment,
        const std::vector<FragmentMetadata> &fragmentList,
        const bool forceWellAnchored) const;

    bool flagDodgyTemplate(FragmentMetadata &orphan, FragmentMetadata &shadow, BamTemplate &bamTemplate) const;
};

inline std::ostream & operator << (std::ostream & os, const TemplateBuilder::BestPairInfo& bestPairInfo)
{
    return os << "BestPairInfo(" <<
        *bestPairInfo.bestPairFragments[0][0] << "-" << *bestPairInfo.bestPairFragments[1][0] << "," <<
        bestPairInfo.resolvedTemplateCount << "rtc, " <<
        bestPairInfo.bestTemplateLogProbability << ":" << bestPairInfo.totalTemplateProbability << "bp:tp, " <<
        bestPairInfo.bestTemplateScore << "bs, " <<
        bestPairInfo.bestPairEditDistance << "bed)";
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_HH
