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
 ** \file TemplateLengthStatistics.hh
 **
 ** \brief Component to encapsulate statistics about template lengths.
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_TEMPLATE_LENGTH_STATISTICS_HH
#define iSAAC_ALIGNMENT_TEMPLATE_LENGTH_STATISTICS_HH

#include <vector>
#include <algorithm>
#include <iostream>
#include <boost/noncopyable.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/constants/constants.hpp>

#include "alignment/FragmentMetadata.hh"
#include "reference/Contig.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Statistics for template lengths
 **
 ** TODO: add support for circular reference
 **
 ** TODO: add support for discordant alignment models
 **/
class TemplateLengthStatistics
{
public:
    enum AlignmentModel
    {
        FFp = 0,
        FRp = 1,
        RFp = 2,
        RRp = 3,
        FFm = 4,
        FRm = 5,
        RFm = 6,
        RRm = 7,
        InvalidAlignmentModel = 8
    };
    TemplateLengthStatistics();

    void clear();


    // Constructor for unit tests
    TemplateLengthStatistics(unsigned min, unsigned max, unsigned median, unsigned lowStdDev, unsigned highStdDev, AlignmentModel m0, AlignmentModel m1,
                             int mateDriftRange, bool stable = true)
        : min_(min)
        , max_(max)
        , median_(-1U)
        , lowStdDev_(lowStdDev)
        , highStdDev_(highStdDev)
        , stable_(stable)
        , mateMin_(-1U)
        , mateMax_(-1U)
    {
        setMedian(median, mateDriftRange);
        bestModels_[0] = m0;
        bestModels_[1] = m1;
    }

    unsigned getMin() const {return min_;}
    unsigned getMax() const {return max_;}
    unsigned getMedian() const {return median_;}
    unsigned getLowStdDev() const {return lowStdDev_;}
    unsigned getHighStdDev() const {return highStdDev_;}

    AlignmentModel getBestModel(unsigned i) const {
        return i > 1 ? InvalidAlignmentModel : static_cast<AlignmentModel>(bestModels_[i]);
    }
    bool isStable() const {return stable_;}
    unsigned getMateMin() const {return mateMin_;}
    unsigned getMateMax() const {return mateMax_;}

    enum CheckModelResult
    {
        Oversized = 0,
        Undersized = 1,
        Nominal = 2,
        NoMatch = 3,
        CheckModelLast
    };

    template <typename FragmenT> CheckModelResult checkModel(
        const FragmenT &f1, const FragmenT &f2) const
    {
        if (f1.getContigId() == f2.getContigId())
        {
            const AlignmentModel model = alignmentModel(f1, f2);
            if (model == bestModels_[0] || model == bestModels_[1])
            {
                const unsigned long length = getLength(f1, f2);
                return (length > max_) ? Oversized : (length < min_) ? Undersized : Nominal;
            }
        }

        return NoMatch;
    }

    /// Determines if the mates should have opposite orientations
    //bool oppositeOrientations(unsigned i) const
    //{
    //    return (bestModels_[i] & 1) ^ ((bestModels_[i] >> 1) & 1);
    //}
    /// Determine the expected orientation of the mate (true is reverse)
    bool mateOrientation(unsigned readIndex, bool reverse) const;
    /// Determine the minimal mate position for the given fragment
    long mateMinPosition(unsigned readIndex, bool reverse, long position, const unsigned *readLengths) const;
    /// Determine the maximal mate position for the given fragment
    long mateMaxPosition(unsigned readIndex, bool reverse, long position, const unsigned *readLengths) const;
    /**
     ** \brief checks if two fragments match the statistical model
     **
     ** Verifies that both the alignment model and the length are as expected
     **/
    bool matchModel(const FragmentMetadata &f1, const FragmentMetadata &f2) const;
    /**
     ** \brief calculate the alignment model for a pair of fragments
     **
     ** The model is based on the orientation ('F' or 'R' respectively for
     ** worward and reverse) of the fragments and on their relative position
     ** ('p' if f1<=f2, 'm' otherwise):
     **
     ** - 0: FFp -> Fp
     ** - 1: FRp -> Rp
     ** - 2: RFp -> Rm
     ** - 3: RRp -> Fm
     ** - 4: FFm -> Fm
     ** - 5: FRm -> Rm
     ** - 6: RFm -> Rp
     ** - 7: RRm -> Fp
     **/
    template <typename FragmenT> static AlignmentModel alignmentModel(const FragmenT &f1, const FragmenT &f2)
    {
        if (f1.getContigId() == f2.getContigId())
        {
            const unsigned positionMask = (f1.getPosition() <= f2.getPosition()) ? 0 : 4;
            const unsigned f1Mask = f1.isReverse() ? 2 : 0;
            const unsigned f2Mask = f2.isReverse() ? 1 : 0;
            return static_cast<AlignmentModel>(positionMask | f1Mask | f2Mask);
        }
        return InvalidAlignmentModel;
    }

    template <typename FragmenT> static unsigned long getLength(const FragmenT &f1, const FragmenT &f2)
    {
        ISAAC_ASSERT_MSG(f1.getContigId() == f2.getContigId(), "Both fragments must align to the same contig");
        if (f1.getPosition() < f2.getPosition())
        {
            return std::max<long>(f2.getPosition() + f2.getObservedLength() - f1.getPosition(), f1.getObservedLength());
        }
        else
        {
            return std::max<long>(f1.getPosition() + f1.getObservedLength() - f2.getPosition(), f2.getObservedLength());
        }
    }

    static const std::string &alignmentModelName(AlignmentModel alignmentModel);

    enum AlignmentClass
    {
        Fp = 0,
        Rp = 1,
        Rm = 2,
        Fm = 3,
        InvalidAlignmentClass
    };
    static const std::string &alignmentClassName(AlignmentClass alignmentClass);
    // check the coherence of the alignment class for the two best models
    bool isCoherent() const
    {
        return (bestModels_[0] != bestModels_[1]) &&
            (alignmentClass(bestModels_[0]) == alignmentClass(bestModels_[1]));
    }

    /**
     ** \brief Calculate the class of the alignmentModel (Fp, Rp, Rm or Fm)
     **/
    static AlignmentClass alignmentClass(AlignmentModel model);

    void setMin(unsigned min, int mateDriftRange) {min_ = min; mateMin_ = -1 == mateDriftRange ? min_ : median_ - mateDriftRange;}
    void setMedian(unsigned median, int mateDriftRange)
    {
        median_ = median;
        mateMin_ = -1 == mateDriftRange ? min_ : median_ - mateDriftRange;
        mateMax_ = -1 == mateDriftRange ? max_ : median_ + mateDriftRange;
    }
    void setMax(unsigned max, int mateDriftRange) {max_ = max; mateMax_ = -1 == mateDriftRange ? max_ : median_ + mateDriftRange;}
    void setLowStdDev(unsigned lowStdDev) {lowStdDev_ = lowStdDev;}
    void setHighStdDev(unsigned highStdDev) {highStdDev_ = highStdDev;}
    void setBestModel(AlignmentModel bestModel, unsigned i) {bestModels_[i] = bestModel;}
    void setStable(bool stable) {stable_ = stable;}
    /*
     * \brief enable serialization
     */
    template <class Archive> friend void serialize(Archive &ar, TemplateLengthStatistics &bm, const unsigned int version);

public:
    /// TEMPLATE_LENGTH_THRESHOLD is the longest template length achievable with the existing chemistry. Change this if technology gets better.
    static const unsigned TEMPLATE_LENGTH_THRESHOLD = 50000;
private:
    static const unsigned READS_MAX = 2;

    unsigned min_;
    unsigned max_;
    unsigned median_;
    unsigned lowStdDev_;
    unsigned highStdDev_;
    AlignmentModel bestModels_[2];
    bool stable_;

    // These two are not part of the length statistics. They are speedup helpers for the components that require
    // estimating mate rescue range many times in a row.
    unsigned mateMin_;
    unsigned mateMax_;

    /// Check that at least one of the bestModels_ has the given orientation for the indicated readIndex.
    bool isValidModel(const bool reverse, const unsigned readIndex) const;
    /// Chech if the fragment with the given orientation and readIndex is supposed to come first
    bool firstFragment(const bool reverse, const unsigned readIndex) const;

};

inline std::ostream &operator<<(std::ostream &os, const TemplateLengthStatistics &t)
{
    return os << "TemplateLengthStatistics("
              << t.getMin() << ":"
              << t.getMedian() << ":"
              << t.getMax() << ":"
              << t.getLowStdDev() << ":"
              << t.getHighStdDev() << ", "
              << t.getBestModel(0) << ":"
              << t.getBestModel(1) << ", "
              << t.isStable()  << ")"
              << std::endl;
}


/**
 ** \brief Accumulates and computes template length statistics
 **
 ** TODO: add support for circular reference
 **
 ** TODO: add support for discordant alignment models
 **/
class TemplateLengthDistribution
{
public:
    TemplateLengthDistribution(int mateDriftRange);

    void reserve(const unsigned reserveClusters);
    void unreserve();
    void clear();
    void reset(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList);

    // Constructor for unit tests
    TemplateLengthDistribution(
        unsigned min, unsigned max, unsigned median, unsigned lowStdDev, unsigned highStdDev,
        TemplateLengthStatistics::AlignmentModel m0, TemplateLengthStatistics::AlignmentModel m1,
        bool stable = true)
        : stats_(min, max, median, lowStdDev, highStdDev, m0, m1, stable)
        , mateDriftRange_(-1)
        , templateCount_(0)
        , uniqueCount_(0)
        , count_(0)
    {
    }

    /**
     ** \brief add a template to the statistical model
     **
     ** Considers only the uniquely aligned fragments where both reads are on
     ** the same contig, where each fragment is completely on the contig and
     ** where the corresponding template would be shorter than a pre-decided
     ** length (see defaultTemplateLengthThreshold).
     **
     ** \param fragments the lists of fragments identified for each read
     **
     ** \return true iff the model is stable
     **/
    bool addTemplate(const std::vector<std::vector<FragmentMetadata> > &fragments);
    /// finalize the model after adding all available templates
    bool finalize();

    int getMateDriftRange() const {return mateDriftRange_;}

    const TemplateLengthStatistics &getStatistics() const {return stats_;}
    bool isStable() const {return stats_.isStable();}
protected:
    void setMin(unsigned min) {stats_.setMin(min, mateDriftRange_);}
    void setMedian(unsigned median){stats_.setMedian(median, mateDriftRange_);}
    void setMax(unsigned max) {stats_.setMax(max, mateDriftRange_);}
private:
    static const unsigned READS_MAX = 2;
    static const unsigned UPDATE_FREQUENCY = 10000;
    static const double STANDARD_DEVIATIONS_MAX = 3.0;
    static const double FRAGMENT_LENGTH_CONFIDENCE_INTERVAL;
    static const double FRAGMENT_LENGTH_CONFIDENCE_INTERVAL_1Z;
    static const double LOWER_PERCENT;
    static const double UPPER_PERCENT;
    static const double LOWER_PERCENT_1Z;
    static const double UPPER_PERCENT_1Z;

    TemplateLengthStatistics stats_;

    // -1 - use min_ and max_ for mateMin_ and mateMax_
    // >= 0 - use median_ +- mateDriftRange_ for mateMin_ and mateMax_
    int mateDriftRange_;
    std::vector<unsigned> lengthList_;

    unsigned templateCount_;
    unsigned uniqueCount_;
    // count of the lengths actually used to populate the histograms
    unsigned count_;
    std::vector<std::vector<unsigned> > histograms_;

    void updateStatistics();
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_TEMPLATE_LENGTH_STATISTICS_HH
