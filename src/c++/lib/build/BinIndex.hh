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
 ** \file BinIndex.hh
 **
 ** Helps sorting and duplicate marking on a single alignment bin.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BIN_INDEX_HH
#define iSAAC_BUILD_BIN_INDEX_HH

#include "UnsortedAlignment.hh"

namespace isaac
{
namespace build
{

struct BinIndex
{
    BinIndex(const UnsortedAlignment *alignmentRecord) :
        key_(alignmentRecord), alignmentRecord_(alignmentRecord){}

    struct Key
    {
        Key(const UnsortedAlignment *alignmentRecord) :
            alignmentPos_(alignmentRecord->alignmentPos_)
        {
            if (alignmentRecord->isSingleton())
            {
                rank_.singletonShadow_.shadowBases_ =
                    alignmentRecord->next()->seq[0] |
                    alignmentRecord->next()->seq[1] << 8 |
                    alignmentRecord->next()->seq[2] << 16 |
                    alignmentRecord->next()->seq[3] << 24;
                rank_.singletonShadow_.singletonSras_ = alignmentRecord->sras_;
                templateLength_ = 0;
            }
            else if (alignmentRecord->isShadow())
            {
                rank_.singletonShadow_.shadowBases_ =
                    alignmentRecord->next()->seq[0] |
                    alignmentRecord->next()->seq[1] << 8 |
                    alignmentRecord->next()->seq[2] << 16 |
                    alignmentRecord->next()->seq[3] << 24;
                rank_.singletonShadow_.singletonSras_ = alignmentRecord->mateSras_;
                templateLength_ = -1;
            }
            else //ohterwise it's either NMNM, chimera, or one of the pair types
            {
                rank_.pair_.srassum_ = alignmentRecord->sras_ + alignmentRecord->mateSras_;
                rank_.pair_.pras_ = alignmentRecord->pras_;
                templateLength_ = alignmentRecord->templateLength_;
            }
        }

        unsigned long alignmentPos_;
        int templateLength_;
        union
        {
            // TODO: see if this can be shortened to 32 bits
            unsigned long value_;
            struct
            {
                // if pras_ is 0 or equals, srassum_ will help rank the pair
                unsigned srassum_;
                // high-order bits in value_ give precedence to pair alignment score
                // over the srassum for pairs that have it
                unsigned pras_;
            } pair_;
            struct
            {
                // 8 first bases of the shadow in a singleton/shadow pair
                // somewhat allow identifying PCR duplicates
                unsigned shadowBases_;
                unsigned singletonSras_;
            } singletonShadow_;
        } rank_;

        bool operator <(const Key &right) const
        {
            return alignmentPos_ < right.alignmentPos_ ||
                      (alignmentPos_ == right.alignmentPos_ &&
                          (templateLength_ < right.templateLength_ ||
                              (templateLength_ == right.templateLength_ && rank_.value_ < right.rank_.value_)));

        }

        bool operator !=(const Key &right) const
        {
            // rank_ is not considered for key equality and only used to bring the best-ranked on top during sort
            return alignmentPos_ != right.alignmentPos_ ||
                templateLength_ != right.templateLength_;
        }

    } key_;

    const UnsortedAlignment *alignmentRecord_;
};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BIN_INDEX_HH
