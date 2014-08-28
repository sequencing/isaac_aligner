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
 ** \file BinSorter.hh
 **
 ** Contains utilities for fragment index duplicate identification and filtering
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_FRAGMENT_INDEX_FILTERING_HH
#define iSAAC_BUILD_FRAGMENT_INDEX_FILTERING_HH

#include "build/BarcodeBamMapping.hh"
#include "build/FragmentIndex.hh"

namespace isaac
{
namespace build
{

/**
 * \brief Order and compares reads to identify duplicates.
 *
 * \param oneLibrarySample if true, the sample index is used instead of lane-barcode. This ensures that PCR
 *        duplicates from different lanes are caught.
 */
template <bool singleLibrarySamples>
struct FDuplicateFilter
{
    const BarcodeBamMapping::BarcodeSampleIndexMap &barcodeSampleIndex_;
    FDuplicateFilter(const BarcodeBamMapping::BarcodeSampleIndexMap &barcodeSampleIndex):
        barcodeSampleIndex_(barcodeSampleIndex){}
    bool less(const PackedFragmentBuffer &fragments,
                     const FStrandFragmentIndex &left,
                     const FStrandFragmentIndex &right) const
    {
        if (left.fStrandPos_ < right.fStrandPos_)
            {return true;}

        if (left.fStrandPos_ == right.fStrandPos_)
        {
            if (left.mate_.anchor_.value_ < right.mate_.anchor_.value_)
                {return true;}

            if (left.mate_.anchor_.value_ == right.mate_.anchor_.value_)
            {
                if (left.mate_.info_.value_ < right.mate_.info_.value_)
                    {return true;}

                if (left.mate_.info_.value_ == right.mate_.info_.value_)
                {
                    const io::FragmentAccessor &leftFragment = fragments.getFragment(left);
                    const io::FragmentAccessor &rightFragment = fragments.getFragment(right);
                    const unsigned long leftLibrary = singleLibrarySamples ?
                        barcodeSampleIndex_.at(leftFragment.barcode_) : leftFragment.barcode_;
                    const unsigned long rightLibrary = singleLibrarySamples ?
                        barcodeSampleIndex_.at(rightFragment.barcode_) : rightFragment.barcode_;
                    //same library must be grouped together as we dupe-remove only within the library
                    if (leftLibrary < rightLibrary)
                        {return true;}

                    if (leftLibrary == rightLibrary)
                    {
                        // we want higher alignment score on top
                        if (left.duplicateClusterRank_ > right.duplicateClusterRank_)
                            {return true;}

                        if (left.duplicateClusterRank_ == right.duplicateClusterRank_)
                        {
                            if (leftFragment.tile_ * INSANELY_HIGH_NUMBER_OF_CLUSTERS_PER_TILE + leftFragment.clusterId_ <
                                rightFragment.tile_ * INSANELY_HIGH_NUMBER_OF_CLUSTERS_PER_TILE + rightFragment.clusterId_)
                                {return true;}
                        }
                    }
                }
            }
        }

        return false;
    }
    bool equal_to(const PackedFragmentBuffer &fragments,
                         const FStrandFragmentIndex &left,
                         const FStrandFragmentIndex &right) const
    {
        if (left.fStrandPos_ == right.fStrandPos_ &&
            left.mate_.anchor_.value_ == right.mate_.anchor_.value_ &&
            left.mate_.info_.value_ == right.mate_.info_.value_)
        {
            const io::FragmentAccessor &leftFragment = fragments.getFragment(left);
            const io::FragmentAccessor &rightFragment = fragments.getFragment(right);

            // In a weird case when both ends of a pair are facing the same way and align at the same position,
            // this condition prevents one being discarded as duplicated of another
            if (leftFragment.tile_ != rightFragment.tile_ ||
                leftFragment.clusterId_ != rightFragment.clusterId_)
            {
                const unsigned long leftLibrary = singleLibrarySamples ?
                    barcodeSampleIndex_.at(leftFragment.barcode_) : leftFragment.barcode_;
                const unsigned long rightLibrary = singleLibrarySamples ?
                    barcodeSampleIndex_.at(rightFragment.barcode_) : rightFragment.barcode_;

                return leftLibrary == rightLibrary;
            }
        }
        return false;
    }
};

template <bool singleLibrarySamples>
struct RSDuplicateFilter
{
    const BarcodeBamMapping::BarcodeSampleIndexMap &barcodeSampleIndex_;
    RSDuplicateFilter(const BarcodeBamMapping::BarcodeSampleIndexMap &barcodeSampleIndex):
        barcodeSampleIndex_(barcodeSampleIndex){}
    bool less(const PackedFragmentBuffer &fragments,
                     const RStrandOrShadowFragmentIndex &left,
                     const RStrandOrShadowFragmentIndex &right) const
    {
        if (left.anchor_.value_ < right.anchor_.value_)
            {return true;}

        if (left.anchor_.value_ == right.anchor_.value_)
        {
            if (left.mate_.anchor_.value_ < right.mate_.anchor_.value_)
                {return true;}

            if (left.mate_.anchor_.value_ == right.mate_.anchor_.value_)
            {
                if (left.mate_.info_.value_ < right.mate_.info_.value_)
                    {return true;}

                if (left.mate_.info_.value_ == right.mate_.info_.value_)
                {
                    const io::FragmentAccessor &leftFragment = fragments.getFragment(left);
                    const io::FragmentAccessor &rightFragment = fragments.getFragment(right);
                    const unsigned long leftLibrary = singleLibrarySamples ?
                        barcodeSampleIndex_.at(leftFragment.barcode_) : leftFragment.barcode_;
                    const unsigned long rightLibrary = singleLibrarySamples ?
                        barcodeSampleIndex_.at(rightFragment.barcode_) : rightFragment.barcode_;
                    //same library must be grouped together as we dupe-remove only within the library
                    if (leftLibrary < rightLibrary)
                        {return true;}

                    if (leftLibrary == rightLibrary)
                    {
                        // we want higher alignment score on top
                        if (left.duplicateClusterRank_ > right.duplicateClusterRank_)
                            {return true;}

                        if (left.duplicateClusterRank_ == right.duplicateClusterRank_)
                        {
                            if (fragments.getFragment(left).tile_ * INSANELY_HIGH_NUMBER_OF_CLUSTERS_PER_TILE + fragments.getFragment(left).clusterId_ <
                                fragments.getFragment(right).tile_ * INSANELY_HIGH_NUMBER_OF_CLUSTERS_PER_TILE + fragments.getFragment(right).clusterId_)
                            {
                                return true;
                            }
                        }
                    }
                }
            }
        }

        return false;
    }
    bool equal_to(const PackedFragmentBuffer &fragments,
                         const RStrandOrShadowFragmentIndex &left,
                         const RStrandOrShadowFragmentIndex &right) const
    {
        if (left.anchor_.value_ == right.anchor_.value_ &&
            left.mate_.anchor_.value_ == right.mate_.anchor_.value_ &&
            left.mate_.info_.value_ == right.mate_.info_.value_)
        {
            const io::FragmentAccessor &leftFragment = fragments.getFragment(left);
            const io::FragmentAccessor &rightFragment = fragments.getFragment(right);

            // In a weird case when both ends of a pair are facing the same way and align at the same position,
            // this condition prevents one being discarded as duplicated of another
            if (leftFragment.tile_ != rightFragment.tile_ ||
                leftFragment.clusterId_ != rightFragment.clusterId_)
            {
                const unsigned long leftLibrary = singleLibrarySamples ?
                    barcodeSampleIndex_.at(leftFragment.barcode_) : leftFragment.barcode_;
                const unsigned long rightLibrary = singleLibrarySamples ?
                    barcodeSampleIndex_.at(rightFragment.barcode_) : rightFragment.barcode_;

                return leftLibrary == rightLibrary;
            }
        }
        else
        {
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(leftFragment.clusterId_, "Not a duplicate " << left << ":" << right);
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(rightFragment.clusterId_, "Not a duplicate " << left << ":" << right);
        }
        return false;
    }

};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_FRAGMENT_INDEX_FILTERING_HH
