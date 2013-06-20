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
 ** TODO: add support for discorant alignment models
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
    TemplateLengthStatistics(int mateDriftRange);

    void reserve(const unsigned reserveClusters);
    void unreserve();
    void clear();
    void reset(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList);

    void setGenome(
        const std::vector<reference::Contig> &contigList,
        const flowcell::ReadMetadataList &readMetadataList);


    // Constructor for unit tests
    TemplateLengthStatistics(unsigned min, unsigned max, unsigned median, unsigned lowStdDev, unsigned highStdDev, AlignmentModel m0, AlignmentModel m1,
                             bool stable = true)
        : mateDriftRange_(-1)
        , min_(min)
        , max_(max)
        , median_(median)
        , lowStdDev_(lowStdDev)
        , highStdDev_(highStdDev)
        , stable_(stable)
        , mateMin_(min_)
        , mateMax_(max_)
        , rogCorrection_(0.0)
    {
        rogCorrectionList_[0] = 0.0;
        rogCorrectionList_[1] = 0.0;
        bestModels_[0] = m0;
        bestModels_[1] = m1;
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
     ** \return true iif the model is stable
     **/
    bool addTemplate(const std::vector<std::vector<FragmentMetadata> > &fragments);
    /// finalize the model after adding all available templates
    bool finalize();
    unsigned getMin() const {return min_;}
    unsigned getMax() const {return max_;}
    unsigned getMedian() const {return median_;}
    unsigned getLowStdDev() const {return lowStdDev_;}
    unsigned getHighStdDev() const {return highStdDev_;}
    unsigned getMateMin() const {return mateMin_;}
    unsigned getMateMax() const {return mateMax_;}

    AlignmentModel getBestModel(unsigned i) const {
        return i > 1 ? InvalidAlignmentModel : static_cast<AlignmentModel>(bestModels_[i]);
    }
    bool isStable() const {return stable_;}
    int getMateDriftRange() const {return mateDriftRange_;}

    enum CheckModelResult
    {
        Oversized = 0,
        Undersized = 1,
        Nominal = 2,
        NoMatch = 3,
        CheckModelLast
    };
    CheckModelResult checkModel(
        const FragmentMetadata &f1, const FragmentMetadata &f2) const;
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
    static AlignmentModel alignmentModel(const FragmentMetadata &f1, const FragmentMetadata &f2);
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

    double getReadRogCorrection(const unsigned readIndex) const
    {
        ISAAC_ASSERT_MSG(sizeof(rogCorrectionList_) > readIndex, "Only up to 2 reads supported");
        return rogCorrectionList_[readIndex];
    }
    double getRogCorrection() const {return rogCorrection_;}
protected:
    void setMin(unsigned min) {min_ = min; mateMin_ = -1 == mateDriftRange_ ? min_ : median_ - mateDriftRange_;}
    void setMedian(unsigned median)
    {
        median_ = median;
        mateMin_ = -1 == mateDriftRange_ ? min_ : median_ - mateDriftRange_;
        mateMax_ = -1 == mateDriftRange_ ? max_ : median_ + mateDriftRange_;
    }
    void setMax(unsigned max) {max_ = max; mateMax_ = -1 == mateDriftRange_ ? max_ : median_ + mateDriftRange_;}
    void setLowStdDev(unsigned lowStdDev) {lowStdDev_ = lowStdDev;}
    void setHighStdDev(unsigned highStdDev) {highStdDev_ = highStdDev;}
    void setBestModel(AlignmentModel bestModel, unsigned i) {bestModels_[i] = bestModel;}
    void setStable(bool stable) {stable_ = stable;}
private:
    static const unsigned readsMax_ = 2;
    static const unsigned defaultTemplateLengthThreshold = 50000;
    static const unsigned updateFrequency = 10000;
    static const double defaultNumberOfStandardDeviations;
    static const double fragmentLengthConfidenceInterval;
    static const double fragmentLengthConfidenceInterval1Z;
    static const double lowerPercent;
    static const double upperPercent;
    static const double lowerPercent1Z;
    static const double upperPercent1Z;
    // -1 - use min_ and max_ for mateMin_ and mateMax_
    // >= 0 - use median_ +- mateDriftRange_ for mateMin_ and mateMax_
    int mateDriftRange_;
    unsigned min_;
    unsigned max_;
    unsigned median_;
    unsigned lowStdDev_;
    unsigned highStdDev_;
    AlignmentModel bestModels_[2];
    std::vector<unsigned> lengthList_;
    bool stable_;
    unsigned mateMin_;
    unsigned mateMax_;

    unsigned templateCount_;
    unsigned uniqueCount_;
    // count of the lengths actually used to populate the histograms
    unsigned count_;
    std::vector<std::vector<unsigned> > histograms_;

    /// Rest-of-genome correction for individual fragments
    double rogCorrectionList_[readsMax_];
    /// Rest-of-genome correction for the template when all fragments match
    double rogCorrection_;

    void updateStatistics();
    unsigned long getLength(const FragmentMetadata &f1, const FragmentMetadata &f2) const;
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
              << t.getMateDriftRange() << ", "
              << t.isStable()  << ")"
              << std::endl;
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_TEMPLATE_LENGTH_STATISTICS_HH
