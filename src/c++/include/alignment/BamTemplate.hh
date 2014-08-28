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
 ** \file BamTemplate.hh
 **
 ** \brief DNA/RNA sequence composed of one or several Fragment(s), as defined
 ** by the SAM Format Specification
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_TEMPLATE_HH
#define iSAAC_ALIGNMENT_TEMPLATE_HH

#include <iostream>
#include <numeric>
#include <vector>

#include "alignment/FragmentMetadata.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Container to encapsulate all the data and metadata associated to a DNA/RNA BamTemplate.
 **
 ** \sa Fragment
 ** \sa FragmentId
 ** \sa Cigar
 ** \sa FragmentMetadata
 **/
class BamTemplate: boost::noncopyable
{
public:
    /**
     ** \brief Constructor using a shared buffer for fragment CIGAR
     **
     ** Note: it is the responsibility of the caller code to ensure the validity
     ** of the cigarBuffer throughout the useful existence of the BamTemplate
     ** instance.
     **/
    BamTemplate(const std::vector<unsigned> &cigarBuffer);
    /**
     ** \brief Initialization for a given cluster
     **
     ** Create the appropriate unaligned FragmentMetadata for each read in the
     ** cluster.
     **/
    void initialize(const flowcell::ReadMetadataList &tileReads, const Cluster &cluster);

    unsigned getMismatchCount() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::getMismatchCount, _2)));
    }

    unsigned getQuality() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0U,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::getQuality, _2)));
    }

    unsigned getEditDistance() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::getEditDistance, _2)));
    }

    unsigned getTotalReadLength() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::getReadLength, _2)));
    }

    bool isUnanchored() const
    {
        return fragmentMetadataList_.end() ==
            std::find_if(fragmentMetadataList_.begin(), fragmentMetadataList_.end(),
                         boost::bind(&FragmentMetadata::getAlignmentScore, _1) != 0);
    }

    unsigned getFragmentCount() const {return fragmentMetadataList_.size();}
    const FragmentMetadata &getFragmentMetadata(unsigned fragmentIndex) const {return fragmentMetadataList_.at(fragmentIndex);}
    const FragmentMetadata &getMateFragmentMetadata(const FragmentMetadata &mate) const {return fragmentMetadataList_.at(getFragmentCount() - 1 - mate.getReadIndex());}
    FragmentMetadata &getFragmentMetadata(unsigned fragmentIndex) {return fragmentMetadataList_[fragmentIndex];}
    FragmentMetadata &getMateFragmentMetadata(const FragmentMetadata &mate) {return fragmentMetadataList_.at(getFragmentCount() - 1 - mate.getReadIndex());}
    const std::vector<unsigned> &getCigarBuffer() const {return cigarBuffer_;}
    unsigned getAlignmentScore() const {return alignmentScore_;}
    bool hasAlignmentScore() const {return -1U != alignmentScore_;}
    void setAlignmentScore(unsigned alignmentScore) {alignmentScore_ = alignmentScore;}
    // TODO: implement passes filter
    bool getPassesFilter() const {return fragmentMetadataList_[0].getCluster().getPf();}
    void setProperPair(const bool properPair) {properPair_ = properPair;}
    bool isProperPair() const {return properPair_;}
    bool filterLowQualityFragments(const unsigned mapqThreshold);
private:
    friend std::ostream &operator<<(std::ostream &os, const BamTemplate &bamTemplate);

    //std::vector<Fragment> fragmentList_;
    std::vector<FragmentMetadata> fragmentMetadataList_;
    const std::vector<unsigned> &cigarBuffer_;

    /**
     ** This depends on the all the pLogCorrect values for all the possible
     ** alignments for this template across the whole reference. It also takes into
     ** account the rest-of-genome correction. Value of -1U indicates unknown alignment score.
     **/
    unsigned alignmentScore_;
    bool properPair_;
};

inline std::ostream &operator<<(std::ostream &os, const BamTemplate &bamTemplate)
{
    return os << "BamTemplate("
         << bamTemplate.fragmentMetadataList_.at(0) << "-" << bamTemplate.fragmentMetadataList_.at(1) << "," <<
         bamTemplate.alignmentScore_ << "as," << bamTemplate.properPair_<< ")";
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_TEMPLATE_HH
