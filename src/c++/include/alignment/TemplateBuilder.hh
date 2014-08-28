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
 ** \file TemplateBuilder.hh
 **
 ** \brief Construction of BamTemplate instances
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_HH
#define iSAAC_ALIGNMENT_TEMPLATE_BUILDER_HH

#include <boost/noncopyable.hpp>

#include "alignment/BamTemplate.hh"
#include "alignment/BandedSmithWaterman.hh"
#include "alignment/Cigar.hh"
#include "alignment/Cluster.hh"
#include "alignment/Match.hh"
#include "alignment/RestOfGenomeCorrection.hh"
#include "alignment/SeedMetadata.hh"
#include "alignment/ShadowAligner.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "flowcell/ReadMetadata.hh"
#include "reference/Contig.hh"

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
	typedef short DodgyAlignmentScore;
	static const DodgyAlignmentScore DODGY_ALIGNMENT_SCORE_UNKNOWN=255;
	static const DodgyAlignmentScore DODGY_ALIGNMENT_SCORE_UNALIGNED=-1;

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
        const bool scatterRepeats,
        const unsigned gappedMismatchesMax,
        const bool avoidSmithWaterman,
        const int gapMatchScore,
        const int gapMismatchScore,
        const int gapOpenScore,
        const int gapExtendScore,
        const int minGapExtendScore,
        const unsigned semialignedGapLimit,
        const DodgyAlignmentScore dodgyAlignmentScore);

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
     ** \return false means template ended up not having a single read aligned anywhere.
     **
     ** This method will initialize the internal template of the builder.
     **
     ** Precondition: the input list of fragments is partitioned by readIndex
     ** and sorted genomic position. This means that the order is (tileId,
     ** clusterId, seedIndex, reverse, contig, position) where 'tileId' and
     ** 'clusterId' are constant.
     **/
    bool buildTemplate(
        const std::vector<reference::Contig> &contigList,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
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
        const RestOfGenomeCorrection &restOfGenomeCorrection,
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
    static const unsigned TRACKED_REPEATS_MAX_ONE_READ = 1000;
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
    class ShadowProbability
    {
        reference::ReferencePosition pos_;
        double logProbability_;
        long observedLength_;
	public:
        ShadowProbability(const FragmentMetadata &shadow) :
            pos_(shadow.getFStrandReferencePosition()), logProbability_(shadow.logProbability),
            observedLength_(shadow.getObservedLength())
        {
            // encode reverse in position to save space.
            // This structure can easily take 5 extra gigabytes if a boolean is introduced to store reverse (SAAC-478)
            pos_.setNeighbors(shadow.isReverse());
        }

        reference::ReferencePosition pos() const {return pos_;}
        double logProbability() const {return logProbability_;}
        long observedLength() const {return observedLength_;}

        bool operator < (const ShadowProbability &that) const
        {
            //if we happen to have fragment and its inversion, let's have the higher probability on top so that it counts towards the total
            return pos_ < that.pos_ ||
                (pos_ == that.pos_ &&
                    (ISAAC_LP_LESS(logProbability_, that.logProbability_) ||
                        (ISAAC_LP_EQUALS(logProbability_, that.logProbability_) &&
                            (observedLength_ < that.observedLength_)
                        )
                    )
                );
        }

        bool operator == (const ShadowProbability &that) const
        {
            return pos_ == that.pos_ && ISAAC_LP_EQUALS(logProbability_, that.logProbability_) &&
                observedLength_ == that.observedLength_;
        }

        bool operator != (const ShadowProbability &that) const { return !(*this == that); }
    };

    // Max number of orphans times the max number of shadows that can be rescued for each orphan plus the max number of shadows that have been discovered during seed matching
    //common::FiniteCapacityVector<ShadowProbability, TRACKED_REPEATS_MAX_ONE_READ * TRACKED_REPEATS_MAX_ONE_READ + TRACKED_REPEATS_MAX_ONE_READ> allShadowProbabilities_[readsMax_];
    std::vector<ShadowProbability> allShadowProbabilities_[readsMax_];

    struct PairProbability
    {
        ShadowProbability r1_;
        ShadowProbability r2_;

        PairProbability(const FragmentMetadata &r1,
                        const FragmentMetadata &r2) : r1_(r1), r2_(r2){}

        bool operator < (const PairProbability &that) const
        {
            // since the sum of log probabilites of reads has to be considered we can't compare reads individually
            return
                (r1_.pos() < that.r1_.pos() || (r1_.pos() == that.r1_.pos() &&
                    (r2_.pos() < that.r2_.pos() || (r2_.pos() == that.r2_.pos() &&
                        //if we happen to have pair and its inversion, let's have the higher probability on top so that it counts towards the total
                        (ISAAC_LP_LESS(that.logProbability(), logProbability()) || (ISAAC_LP_EQUALS(logProbability(), that.logProbability()) &&
                            (r1_.observedLength() < that.r1_.observedLength() || (r1_.observedLength() == that.r1_.observedLength() &&
                                r2_.observedLength() < that.r2_.observedLength()))))))));
        }

        bool operator == (const PairProbability &that) const
        {
            // since the sum of log probabilites of reads has to be considered we can't compare reads individually
            return r1_.pos() == that.r1_.pos() && r2_.pos() == that.r2_.pos() && ISAAC_LP_EQUALS(logProbability(), that.logProbability()) &&
                r1_.observedLength() == that.r1_.observedLength() && r2_.observedLength() == that.r2_.observedLength();
        }

        bool operator != (const PairProbability &that) const {return !(*this == that);}

        double logProbability() const {return r1_.logProbability() + r2_.logProbability();}
    };

    // Max number of orphans times the max number of shadows that can be rescued for each orphan times the number of reads
//    common::FiniteCapacityVector<PairProbability, TRACKED_REPEATS_MAX_ONE_READ * TRACKED_REPEATS_MAX_ONE_READ * readsMax_> allPairProbabilities_;
    std::vector<PairProbability> allPairProbabilities_;



    common::FiniteCapacityVector<FragmentMetadata, TRACKED_REPEATS_MAX_ONE_READ> bestOrphanShadows_[readsMax_];

    typedef std::vector<FragmentMetadata>::const_iterator FragmentIterator;

    struct BestPairInfo
    {
        BestPairInfo() :
            bestTemplateLogProbability(-std::numeric_limits<double>::max()),
            bestTemplateScore(-1),
            resolvedTemplateCount(0), bestPairEditDistance(0),
            totalTemplateProbability(0.0)
        {
        }

        void clear()
        {
            bestTemplateLogProbability = -std::numeric_limits<double>::max();
            bestTemplateScore = -1;
            resolvedTemplateCount = 0;
            bestPairEditDistance = 0;
            totalTemplateProbability = 0.0;

            bestPairFragments[0].clear();
            bestPairFragments[1].clear();
        }

        void init(FragmentIterator bestR1Fragment, FragmentIterator bestR2Fragment)
        {
            clear();
            bestPairFragments[0].push_back(bestR1Fragment);
            bestPairFragments[1].push_back(bestR2Fragment);
        }

        long getBestTemplateLength() const
        {
            if (!resolvedTemplateCount)
            {
                return 0;
            }
            reference::ReferencePosition templateStart = std::min(bestPairFragments[0][0]->getFStrandReferencePosition(),
                                                                  bestPairFragments[1][0]->getFStrandReferencePosition());
            reference::ReferencePosition templateEnd = std::max(bestPairFragments[0][0]->getRStrandReferencePosition(),
                                                                  bestPairFragments[1][0]->getRStrandReferencePosition());
            return templateEnd - templateStart;
        }

        // potentially this needs to hold all the combinations of individual fragment alignments
        typedef isaac::common::FiniteCapacityVector<FragmentIterator, TRACKED_REPEATS_MAX_ONE_READ * TRACKED_REPEATS_MAX_ONE_READ> FragmentIteratorVector;
        FragmentIteratorVector bestPairFragments[readsMax_];
        double bestTemplateLogProbability;
        unsigned long bestTemplateScore;
        unsigned resolvedTemplateCount;
        unsigned bestPairEditDistance;
        double totalTemplateProbability;
    };
    friend std::ostream & operator << (std::ostream & os, const BestPairInfo& bestPairInfo);

    /// Holds the information about pairs obtained via combining alignments from MatchFinder
    BestPairInfo bestCombinationPairInfo_;
    /// Holds the information about the pairs rescued via rescueShadow or buildDisjoinedTemplate
    BestPairInfo bestRescuedPair_;

    /// Helper method to select the best fragment for single-ended runs
    bool pickBestFragment(
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const TemplateLengthStatistics &templateLengthStatistics,
        const std::vector<FragmentMetadata> &fragmentList);
    bool rescueShadow(
        const std::vector<reference::Contig> &contigList,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const TemplateLengthStatistics &templateLengthStatistics);
    /// Helper method to find the best pair of fragments for paired-end runs
    bool pickBestPair(
        const std::vector<reference::Contig> &contigList,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const TemplateLengthStatistics &templateLengthStatistics);

    /// Helper method to locate the best pair of fragments
    void locateBestPair(
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const TemplateLengthStatistics &templateLengthStatistics,
        TemplateBuilder::BestPairInfo &ret) const;
    /// Helper method to build a paired-end template
    bool buildPairedEndTemplate(
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const TemplateLengthStatistics &templateLengthStatistics,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        BestPairInfo &bestPairInfo);
    /// Helper method to produce a template when no pair of fragments match
    bool buildDisjoinedTemplate(
        const std::vector<reference::Contig> &contigList,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const TemplateLengthStatistics &templateLengthStatistics,
        const BestPairInfo &knownBestPair);

    bool scoreDisjoinedTemplate(
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
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
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const TemplateLengthStatistics &templateLengthStatistics,
        const FragmentIterator listFragment,
        const std::vector<FragmentMetadata> &fragmentList,
        const bool forceWellAnchored) const;

    bool flagDodgyTemplate(FragmentMetadata &orphan, FragmentMetadata &shadow, BamTemplate &bamTemplate) const;
    bool flagDodgyTemplate(FragmentMetadata &orphan, BamTemplate &bamTemplate) const;
    static double sumUniqueShadowProbabilities(std::vector<ShadowProbability>::iterator begin, std::vector<ShadowProbability>::iterator end);
    static double sumUniquePairProbabilities(std::vector<PairProbability>::iterator begin, std::vector<PairProbability>::iterator end);

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
