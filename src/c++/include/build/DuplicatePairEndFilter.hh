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
 ** \file DuplicatePairEndFilter.hh
 **
 ** Template for general approach to filtering ends of duplicate pairs.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_DUPLICATE_PAIR_END_FILTER_HH
#define iSAAC_BUILD_DUPLICATE_PAIR_END_FILTER_HH

#include "build/BuildStats.hh"
#include "build/PackedFragmentBuffer.hh"
#include "common/Debug.hh"
#include "io/FragmentIndex.hh"

namespace isaac
{
namespace build
{

template <typename InputT> struct DuplicateFilterTraits
{
//    static bool less(const InputT &left, const InputT &right) {
//        using std::less;
//        return less(left, right);
//    }
//    static bool equal_to(const InputT &left, const InputT &right) {
//        using std::equal_to;
//        return equal_to(left, right);
//    }
    static bool less(const PackedFragmentBuffer &fragments, const InputT &left, const InputT &right) ;
    static bool equal_to(const PackedFragmentBuffer &fragments, const InputT &left, const InputT &right);
    static bool is_reverse();
};

/**
 *
 * \brief This class implements the generic duplicate filtering flow:
 *  1. sort according to the duplicate ranking
 *  2. skip the ones that are duplicates
 *  3. sort the results according to output storage order requirements
 **/
class DuplicatePairEndFilter
{
public:
    DuplicatePairEndFilter(const bool keepDuplicates) : keepDuplicates_(keepDuplicates){}
    template <typename InputIteratorT, typename InsertIteratorT>
    void filterInput(
        PackedFragmentBuffer &fragments,
        InputIteratorT duplicatesBegin,
        InputIteratorT duplicatesEnd,
        BuildStats &buildStats,
        const unsigned binIndex,
        InsertIteratorT results)
    {
        if (duplicatesBegin != duplicatesEnd)
        {
            // reorder them according to duplicate removal rules
            ISAAC_THREAD_CERR << "Sorting duplicates" << std::endl;
            const clock_t startSort = clock();

            std::sort(duplicatesBegin, duplicatesEnd,
                      boost::bind(&DuplicateFilterTraits<typename InputIteratorT::value_type>::less,
                                  boost::ref(fragments), _1, _2));

            ISAAC_THREAD_CERR << "Sorting duplicates" << " done in " << (clock() - startSort) / 1000 << "ms" << std::endl;

            // populate self with the unique fragments
            ISAAC_THREAD_CERR << "Filtering duplicates" << std::endl;
            const clock_t startFilter = clock();


            // Range is guaranteed to be not empty
            unsigned long unique = 0;
            results++ = PackedFragmentBuffer::Index(*duplicatesBegin, fragments.getFragment(*duplicatesBegin));
            ISAAC_THREAD_CERR_DEV_TRACE((boost::format("Selected as the first duplicate best: %s:%s") % *duplicatesBegin %
                fragments.getFragment(*duplicatesBegin)).str());
            for (InputIteratorT it(duplicatesBegin + 1), itLast(duplicatesBegin); duplicatesEnd != it; ++it)
            {
                io::FragmentAccessor &fragment = fragments.getFragment(*it);
                ISAAC_DEV_TRACE_BLOCK(const io::FragmentAccessor &lastFragment = fragments.getFragment(*itLast);)

                if (!DuplicateFilterTraits<typename InputIteratorT::value_type>::equal_to(fragments, *itLast, *it))
                {
                    ISAAC_THREAD_CERR_DEV_TRACE("Selected as a duplicate best:         " << *it << ":" << fragment);
                    ISAAC_THREAD_CERR_DEV_TRACE("Selected as a duplicate best prev:    " << *itLast << ":" << fragments.getFragment(*itLast));
                    results++ = PackedFragmentBuffer::Index(*it, fragment);
                    unique++;
                    itLast = it;
                    buildStats.incrementUniqueFragments(binIndex, fragment.barcode_);
                }
                else if (keepDuplicates_)
                {
                    fragment.flags_.duplicate_ = true;
                    results++ = PackedFragmentBuffer::Index(*it, fragment);
                    ISAAC_THREAD_CERR_DEV_TRACE("Marked as a duplicate:                " << *it << ":" << fragment);
                }
                else
                {
                    ISAAC_THREAD_CERR_DEV_TRACE("Discarded as a duplicate:             " << *it << ":" << fragment);
                }
                buildStats.incrementTotalFragments(binIndex, fragments.getFragment(*it).barcode_);
            }

            ISAAC_THREAD_CERR << "Filtering duplicates"
                << " done in " << (clock() - startFilter) / 1000 << "ms. found " << unique
                << " unique out of " << duplicatesEnd - duplicatesBegin << " fragments" << std::endl;
        }
    }
private:
    const bool keepDuplicates_;
};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_DUPLICATE_PAIR_END_FILTER_HH
