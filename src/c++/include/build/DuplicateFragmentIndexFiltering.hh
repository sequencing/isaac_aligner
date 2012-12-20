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
 ** \file BinSorter.hh
 **
 ** Contains utilities for fragment index duplicate identification and filtering
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_FRAGMENT_INDEX_FILTERING_HH
#define iSAAC_BUILD_FRAGMENT_INDEX_FILTERING_HH

#include "io/FragmentIndex.hh"

namespace isaac
{
namespace build
{

static const unsigned long INSANELY_HIGH_NUMBER_OF_CLUSTERS_PER_TILE = 1000000000;

template <>
struct DuplicateFilterTraits<io::FStrandFragmentIndex>
{
    static bool less(const PackedFragmentBuffer &fragments,
                     const io::FStrandFragmentIndex &left,
                     const io::FStrandFragmentIndex &right)
    {
        if (left.fStrandPos_ < right.fStrandPos_ ||
                (left.fStrandPos_ == right.fStrandPos_ &&
                    (left.mate_.anchor_.value_ < right.mate_.anchor_.value_ ||
                        (left.mate_.anchor_.value_ == right.mate_.anchor_.value_ &&
                            (left.mate_.info_.value_ < right.mate_.info_.value_ ||
                                (left.mate_.info_.value_ == right.mate_.info_.value_ &&
                                    (fragments.getFragment(left).barcode_ < fragments.getFragment(right).barcode_ || //same barcode must be grouped together as we dupe-remove only within the barcode
                                        (fragments.getFragment(left).barcode_ == fragments.getFragment(right).barcode_ &&
                                            (left.duplicateClusterRank_ > right.duplicateClusterRank_ || // we want higher alignment score on top
                                                (left.duplicateClusterRank_ == right.duplicateClusterRank_ &&
                                                    (fragments.getFragment(left).tile_ * INSANELY_HIGH_NUMBER_OF_CLUSTERS_PER_TILE + fragments.getFragment(left).clusterId_ <
                                                        fragments.getFragment(right).tile_ * INSANELY_HIGH_NUMBER_OF_CLUSTERS_PER_TILE + fragments.getFragment(right).clusterId_
                                                    )
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
        )
        {
            return true;
        }

//        if (left.dataOffset_ == right.dataOffset_)
//        {
//            ISAAC_THREAD_CERR << "compare================================================== " <<
//                &left << " against " << &right << std::endl;
//            ISAAC_THREAD_CERR << "left:  " << left << std::endl;
//            ISAAC_THREAD_CERR << "right: " << right << std::endl;
//        }

//gcc 4.3 annd 4.4.5 std::sort seems to create temporary copy and compares it to the existing items which causes this assert to fail.
//        ISAAC_ASSERT_MSG(left.dataOffset_ != right.dataOffset_,
//                         "binary copy of f-strand fragment index encountered. dataOffset_ must be unique!");
        return false;
    }
    static bool equal_to(const PackedFragmentBuffer &fragments,
                         const io::FStrandFragmentIndex &left,
                         const io::FStrandFragmentIndex &right)
    {
        if (left.fStrandPos_ == right.fStrandPos_ &&
            left.mate_.anchor_.value_ == right.mate_.anchor_.value_ &&
            left.mate_.info_.value_ == right.mate_.info_.value_)
        {
            const io::FragmentAccessor &leftFragment = fragments.getFragment(left);
            const io::FragmentAccessor &rightFragment = fragments.getFragment(right);
            return (leftFragment.barcode_ == rightFragment.barcode_ &&
                // In a weird case when both ends of a pair are facing the same way and align at the same position,
                // this condition prevents one being discarded as duplicated of another
                leftFragment.clusterId_ != rightFragment.clusterId_);
        }
        return false;
    }

    static bool is_reverse(){return false;}

};


template <>
struct DuplicateFilterTraits<io::RStrandOrShadowFragmentIndex>
{
    static bool less(const PackedFragmentBuffer &fragments,
                     const io::RStrandOrShadowFragmentIndex &left,
                     const io::RStrandOrShadowFragmentIndex &right)
    {
        if (left.anchor_.value_ < right.anchor_.value_ ||
            (left.anchor_.value_ == right.anchor_.value_ &&
                (left.mate_.anchor_.value_ < right.mate_.anchor_.value_ ||
                    (left.mate_.anchor_.value_ == right.mate_.anchor_.value_ &&
                        (left.mate_.info_.value_ < right.mate_.info_.value_ ||
                            (left.mate_.info_.value_ == right.mate_.info_.value_ &&
                                (fragments.getFragment(left).barcode_ < fragments.getFragment(right).barcode_ || //same barcode must be groupped together as we dupe-remove only within the barcode
                                    (fragments.getFragment(left).barcode_ == fragments.getFragment(right).barcode_ &&
                                        (left.duplicateClusterRank_ > right.duplicateClusterRank_ || // we want higher alignment score on top
                                            (left.duplicateClusterRank_ == right.duplicateClusterRank_ &&
                                                (fragments.getFragment(left).tile_ * INSANELY_HIGH_NUMBER_OF_CLUSTERS_PER_TILE + fragments.getFragment(left).clusterId_ <
                                                    fragments.getFragment(right).tile_ * INSANELY_HIGH_NUMBER_OF_CLUSTERS_PER_TILE + fragments.getFragment(right).clusterId_
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
        {
            return true;
        }

//        if (left.dataOffset_ == right.dataOffset_)
//        {
//            ISAAC_THREAD_CERR << "compare================================================== " <<
//                &left << " against " << &right << std::endl;
//            ISAAC_THREAD_CERR << "left:  " << left << std::endl;
//            ISAAC_THREAD_CERR << "right: " << right << std::endl;
//        }
//gcc 4.3 annd 4.4.5 std::sort seems to create temporary copy and compares it to the existing items which causes this assert to fail.
//        ISAAC_ASSERT_MSG(left.dataOffset_ != right.dataOffset_
//                         ,"binary copy of r-strand fragment index encountered. dataOffset_ must be unique!");
        return false;
    }
    static bool equal_to(const PackedFragmentBuffer &fragments,
                         const io::RStrandOrShadowFragmentIndex &left,
                         const io::RStrandOrShadowFragmentIndex &right)
    {
        if (left.anchor_.value_ == right.anchor_.value_ &&
            left.mate_.anchor_.value_ == right.mate_.anchor_.value_ &&
            left.mate_.info_.value_ == right.mate_.info_.value_)
        {
            const io::FragmentAccessor &leftFragment = fragments.getFragment(left);
            const io::FragmentAccessor &rightFragment = fragments.getFragment(right);
            return (leftFragment.barcode_ == rightFragment.barcode_ &&
                // In a weird case when both ends of a pair are facing the same way and align at the same position,
                // this condition prevents one being discarded as duplicated of another
                leftFragment.clusterId_ != rightFragment.clusterId_);
        }
        return false;
    }

    static bool is_reverse(){return true;}

};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_FRAGMENT_INDEX_FILTERING_HH
