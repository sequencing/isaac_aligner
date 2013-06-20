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
 ** \file TemplateLengthStatistics.cpp
 **
 ** \brief See TemplateLengthStatistics.hh
 ** 
 ** \author Come Raczy
 **/


#include "alignment/Quality.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "alignment/Cigar.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

const double TemplateLengthStatistics::defaultNumberOfStandardDeviations = 3.0;
const double TemplateLengthStatistics::fragmentLengthConfidenceInterval = 
    boost::math::erf(defaultNumberOfStandardDeviations/boost::math::constants::root_two<double>());
const double TemplateLengthStatistics::fragmentLengthConfidenceInterval1Z = boost::math::erf(1.0/boost::math::constants::root_two<double>());
const double TemplateLengthStatistics::lowerPercent = (1.0 - fragmentLengthConfidenceInterval) / 2.0;
const double TemplateLengthStatistics::upperPercent = (1.0 + fragmentLengthConfidenceInterval) / 2.0;
const double TemplateLengthStatistics::lowerPercent1Z = (1.0 - fragmentLengthConfidenceInterval1Z) / 2.0;
const double TemplateLengthStatistics::upperPercent1Z = (1.0 + fragmentLengthConfidenceInterval1Z) / 2.0;

TemplateLengthStatistics::TemplateLengthStatistics(int mateDriftRange)
    : mateDriftRange_(mateDriftRange)
    , min_(0)
    , max_(0)
    , median_(0)
    , lowStdDev_(0)
    , highStdDev_(0)
    , stable_(false)
    , templateCount_(0)
    , uniqueCount_(0)
    , count_(0)
    , histograms_(InvalidAlignmentModel)
    , rogCorrection_(0.0)
{
    bestModels_[0] = InvalidAlignmentModel;
    bestModels_[1] = InvalidAlignmentModel;
    std::fill(rogCorrectionList_, rogCorrectionList_ + sizeof(rogCorrectionList_) / sizeof(rogCorrectionList_[0]), 0.0);

}
void TemplateLengthStatistics::reserve(const unsigned reserveClusters)
{
    ISAAC_ASSERT_MSG(reserveClusters, "reserveClusters must not be 0");
    std::for_each(histograms_.begin(), histograms_.end(), boost::bind(&std::vector<unsigned>::reserve, _1, reserveClusters));
    lengthList_.reserve(reserveClusters * 2);
}

void TemplateLengthStatistics::unreserve()
{
    std::vector<std::vector<unsigned> >().swap(histograms_);
    std::vector<unsigned>().swap(lengthList_);
}

inline void getRestOfGenomeCorrection(
    const std::vector<reference::Contig> &contigList,
    const flowcell::ReadMetadataList &readMetadataList,
    double rogCorrectionList[2])
{
    const size_t genomeLength = reference::genomeLength(contigList);
    BOOST_FOREACH(const flowcell::ReadMetadata &readMetadata, readMetadataList)
    {
        rogCorrectionList[readMetadata.getIndex()] = Quality::restOfGenomeCorrection(genomeLength, readMetadata.getLength());
        // can't have 0.0 in it as it turns the alignment score into 0
        rogCorrectionList[readMetadata.getIndex()] = std::max(rogCorrectionList[readMetadata.getIndex()], std::numeric_limits<double>::min());
    }
}

void TemplateLengthStatistics::clear()
{
    min_ = 0;
    max_ = 0;
    median_ = 0;
    lowStdDev_ = 0;
    highStdDev_ = 0;
    stable_ = false;
    templateCount_ = 0;
    uniqueCount_ = 0;
    count_ = 0;
    bestModels_[0] = InvalidAlignmentModel;
    bestModels_[1] = InvalidAlignmentModel;
    std::for_each(histograms_.begin(), histograms_.end(), std::mem_fun_ref(&std::vector<unsigned>::clear));
    lengthList_.clear();
}

void TemplateLengthStatistics::reset(const std::vector<reference::Contig> &contigList,
                                     const flowcell::ReadMetadataList &readMetadataList)
{
    clear();
    setGenome(contigList, readMetadataList);
}

void TemplateLengthStatistics::setGenome(const std::vector<reference::Contig> &contigList,
                                         const flowcell::ReadMetadataList &readMetadataList)
{
    getRestOfGenomeCorrection(contigList, readMetadataList, rogCorrectionList_);
    rogCorrection_ = Quality::restOfGenomeCorrection(reference::genomeLength(contigList),
                                                     flowcell::getTotalReadLength(readMetadataList));
    // can't have 0.0 in it as it turns the alignment score into 0
    rogCorrection_ = std::max(rogCorrection_, std::numeric_limits<double>::min());

}
bool TemplateLengthStatistics::addTemplate(const std::vector<std::vector<FragmentMetadata> > &fragments)
{
    ISAAC_ASSERT_MSG(2 == fragments.size(), "Maximum of two fragments per template is supported");
    // discard templates where at least on fragment didn't align
    if (fragments[0].empty() || fragments[1].empty())
    {
        return isStable();
    }
    ++templateCount_;
    // discard templates where the alignment is not unique on both fragments
    if ((1 < fragments[0].size()) || (1 < fragments[1].size()))
    {
        return isStable();
    }
    ++uniqueCount_;
    // discard templates that span across several contigs
    if (fragments[0][0].contigId != fragments[1][0].contigId)
    {
        return isStable();
    }
    // discard fragments that ar not completely contained inside the contig
    // this is identified by the presence of leading or trailing inserts 
    assert(fragments[0][0].cigarBuffer == fragments[1][0].cigarBuffer);
    const std::vector<unsigned> &cigarBuffer = *fragments[0][0].cigarBuffer;
    for (size_t i = 0; 2 > i; ++i)
    {
        const unsigned firstOp = cigarBuffer[fragments[i][0].cigarOffset];
        const unsigned lastOp = cigarBuffer[fragments[i][0].cigarOffset + fragments[i][0].cigarLength - 1];
        if ((firstOp & 0xF) == Cigar::INSERT || (lastOp & 0xF) == Cigar::INSERT)
        {
            return isStable();
        }
    }
    // calculate the lenght of the template
    const unsigned long length = getLength(fragments[0][0], fragments[1][0]);
    // discard excessingly long templates
    if (length > defaultTemplateLengthThreshold)
    {
        return isStable();
    }
    // update the histogram for the appropriate alignement model
    const AlignmentModel am = alignmentModel(fragments[0][0], fragments[1][0]);
    if (InvalidAlignmentModel != am)
    {
        histograms_[am].push_back(length);
        // calculate the alignment statistics if appropriate
        ++count_;
        if (0 == (count_ % updateFrequency))
        {
            // save the old statistics
            const unsigned oldMin = min_;
            const unsigned oldMax = max_;
            const unsigned oldMedian = median_;
            const unsigned oldLowStdDev = lowStdDev_;
            const unsigned oldHighStdDev = highStdDev_;
            // get the new statistics
            updateStatistics();
            // check if we reached stability
            if (oldMin == min_ &&
                oldMedian == median_ &&
                oldMax == max_ &&
                oldLowStdDev == lowStdDev_ &&
                oldHighStdDev == highStdDev_)
            {
                stable_ = true;
            }
        }
    }
    return isStable();
}

bool TemplateLengthStatistics::finalize()
{
    // save the old statistics
    const unsigned oldMin = min_;
    const unsigned oldMax = max_;
    const unsigned oldMedian = median_;
    const unsigned oldLowStdDev = lowStdDev_;
    const unsigned oldHighStdDev = highStdDev_;
    // get the new statistics
    updateStatistics();
    // check if we reached stability
    if (oldMin == min_ &&
        oldMedian == median_ &&
        oldMax == max_ &&
        oldLowStdDev == lowStdDev_ &&
        oldHighStdDev == highStdDev_)
    {
        stable_ = true;
    }
    return isStable();
}

unsigned long TemplateLengthStatistics::getLength(const FragmentMetadata &f1, const FragmentMetadata &f2) const
{
    assert(f1.contigId == f2.contigId);
    if (f1.position < f2.position)
    {
        return std::max(f2.observedLength + f2.position - f1.position, f1.observedLength);
    }
    else
    {
        return std::max(f1.observedLength + f1.position - f2.position, f2.observedLength);
    }
}

TemplateLengthStatistics::CheckModelResult TemplateLengthStatistics::checkModel(const FragmentMetadata &f1, const FragmentMetadata &f2) const
{
    if (f1.contigId == f2.contigId)
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

bool TemplateLengthStatistics::matchModel(const FragmentMetadata &f1, const FragmentMetadata &f2) const
{

    const unsigned long length = getLength(f1, f2);
    const AlignmentModel model = alignmentModel(f1, f2);
    const bool ret = (length <= max_ + defaultTemplateLengthThreshold) && ((model == bestModels_[0]) || (model == bestModels_[1]));
    ISAAC_THREAD_CERR_DEV_TRACE((boost::format("TemplateLengthStatistics::matchModel %d: %s-%s (len: %u)") % ret %
                                 f1 % f2 % length).str());
    return ret;

}

TemplateLengthStatistics::AlignmentModel TemplateLengthStatistics::alignmentModel(const FragmentMetadata &f1, const FragmentMetadata &f2)
{
    if (f1.contigId == f2.contigId)
    {
        const unsigned positionMask = (f1.position <= f2.position) ? 0 : 4;
        const unsigned f1Mask = f1.reverse ? 2 : 0;
        const unsigned f2Mask = f2.reverse ? 1 : 0;
        return static_cast<AlignmentModel>(positionMask | f1Mask | f2Mask);
    }
    return InvalidAlignmentModel;
}

TemplateLengthStatistics::AlignmentClass TemplateLengthStatistics::alignmentClass(const AlignmentModel model)
{
    assert(model < 8);
    return static_cast<AlignmentClass>((model < 4) ? model : ((~model) & 3));
}

const static std::vector<std::string> classNames = boost::assign::list_of("F+")("R+")("R-")("F-")("unknown");
const std::string &TemplateLengthStatistics::alignmentClassName(const AlignmentClass alignmentClass)
{
    return (classNames.size() > size_t(alignmentClass)) ? classNames.at(alignmentClass) : classNames.back();
}

const static std::vector<std::string> modelNames = boost::assign::list_of
        ("FF+")("FR+")("RF+")("RR+")("FF-")("FR-")("RF-")("RR-")("unknown");
const std::string &TemplateLengthStatistics::alignmentModelName(const AlignmentModel alignmentModel)
{
    return (modelNames.size() > size_t(alignmentModel)) ? modelNames.at(alignmentModel) : modelNames.back();
}

void TemplateLengthStatistics::updateStatistics()
{
    // save the old statistics
    const unsigned oldMin = min_;
    const unsigned oldMax = max_;
    const unsigned oldMedian = median_;
    const unsigned oldLowStdDev = lowStdDev_;
    const unsigned oldHighStdDev = highStdDev_;
    const AlignmentModel oldBestModels[] = {bestModels_[0], bestModels_[1]};
    // find the two best alignment models
    bestModels_[0] = histograms_[1].size() <= histograms_[0].size() ? FFp : FRp;
    bestModels_[1] = static_cast<AlignmentModel>((bestModels_[0] + 1) % 2);
    for (size_t i = 2; histograms_.size() > i; ++i)
    {
        if (histograms_[i].size() > histograms_[bestModels_[0]].size())
        {
            bestModels_[1] = bestModels_[0];
            bestModels_[0] = static_cast<AlignmentModel>(i);
        }
        else if (histograms_[i].size() > histograms_[bestModels_[1]].size())
        {
            bestModels_[1] = static_cast<AlignmentModel>(i);
        }
    }

    // TODO: report incoherent models
    if (alignmentClass(bestModels_[0]) != alignmentClass(bestModels_[1]))
    {
        ISAAC_THREAD_CERR << "Incoherent alignment models:"
            << bestModels_[0] << " (" << alignmentModelName(bestModels_[0]) << "->"
                << alignmentClassName(alignmentClass(bestModels_[0])) << "):"
            << bestModels_[1] << " (" << alignmentModelName(bestModels_[1]) << "->"
                << alignmentClassName(alignmentClass(bestModels_[1])) << ")" << std::endl;
    }
    // sort all the lengths found for the two best models
    lengthList_.clear();
    lengthList_.insert(lengthList_.end(), histograms_[bestModels_[0]].begin(), histograms_[bestModels_[0]].end());
    lengthList_.insert(lengthList_.end(), histograms_[bestModels_[1]].begin(), histograms_[bestModels_[1]].end());
    std::sort(lengthList_.begin(), lengthList_.end());
    // compute the statistics
    setMin(lengthList_.empty() ? 0 : lengthList_[unsigned(lengthList_.size() * lowerPercent)]);
    setMedian(lengthList_.empty() ? defaultTemplateLengthThreshold/2 : lengthList_[unsigned(lengthList_.size() * 0.5)]);
    setMax(lengthList_.empty() ? defaultTemplateLengthThreshold : lengthList_[unsigned(lengthList_.size() * upperPercent)]);
    lowStdDev_ =  lengthList_.empty() ? median_ : (median_ - lengthList_[unsigned(lengthList_.size() * lowerPercent1Z)]);
    highStdDev_ = lengthList_.empty() ? median_ : (lengthList_[unsigned(lengthList_.size() * upperPercent1Z)] - median_);

    // this has to be cleared when unused. Otherwise, copying template length statistics causes unwanted memory allocations
    lengthList_.clear();
    // check if we reached stability
    if (oldMin == min_ &&
        oldMedian == median_ &&
        oldMax == max_ &&
        oldLowStdDev == lowStdDev_ &&
        oldHighStdDev == highStdDev_ &&
        oldBestModels[0] == bestModels_[0] &&
        oldBestModels[1] == bestModels_[1])
    {
        stable_ = true;
    }
}

/// Check that one of the two best models has the given orientation for the readIndex
bool TemplateLengthStatistics::isValidModel(const bool reverse, const unsigned readIndex) const
{
    ISAAC_ASSERT_MSG(readsMax_ > readIndex, "Invalid read index");
    // bit 0 is for read 2; bit 1 is for read 1
    const unsigned shift = (readIndex + 1) % 2;
    return (reverse == ((bestModels_[0] >> shift) & 1)) || (reverse == ((bestModels_[1] >> shift) & 1));
}

bool TemplateLengthStatistics::firstFragment(const bool reverse, const unsigned readIndex) const
{
    ISAAC_ASSERT_MSG(readsMax_ > readIndex, "Invalid read index");
    // bit 0 is for read 2; bit 1 is for read 1
    const unsigned shift = (readIndex + 1) % 2;
    for (unsigned i = 0; 2 > i; ++i)
    {
        if (reverse == ((bestModels_[i] >> shift) & 1))
        {
            return ((bestModels_[i] >> 2) & 1)== readIndex;
        }
    }
    assert(false); // shouldn't get there
}

bool TemplateLengthStatistics::mateOrientation(unsigned readIndex, bool reverse) const
{
    // bit 0 is for readIndex 1 and vice-versa
    const unsigned shift = (readIndex + 1) % 2;
    for (unsigned i = 0; 2 > i; ++i)
    {
        if (reverse == ((bestModels_[i] >> shift) & 1))
        {
            return (bestModels_[i] >> readIndex) & 1;
        }
    }
    // for discorant models, return the expected orientation in the best model
    return (bestModels_[0] >> readIndex) & 1;
}

long TemplateLengthStatistics::mateMinPosition(
    const unsigned readIndex,
    const bool reverse,
    const long position,
    const unsigned *readLengths) const
{
    if (!isValidModel(reverse, readIndex))
    {
        return position;
    }
    if (firstFragment(reverse, readIndex))
    {
        return position + mateMin_ - readLengths[(readIndex + 1) % 2];
    }
    else
    {
        return position - mateMax_ + readLengths[readIndex];
    }
}

long TemplateLengthStatistics::mateMaxPosition(
    const unsigned readIndex,
    const bool reverse,
    const long position,
    const unsigned *readLengths) const
{
    if (!isValidModel(reverse, readIndex))
    {
        return position;
    }
    if (firstFragment(reverse, readIndex))
    {
        return position + mateMax_ - readLengths[(readIndex + 1) % 2];
    }
    else
    {
        return position - mateMin_ + readLengths[readIndex];
    }
}

} // namespace alignment
} // namespace isaac
