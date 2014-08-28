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
 ** \file NotAFilter.hh
 **
 ** Template for general approach to filtering ends of duplicate pairs.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_NOT_A_FILTER_HH
#define iSAAC_BUILD_NOT_A_FILTER_HH

#include "build/BuildStats.hh"
#include "build/FragmentIndex.hh"
#include "build/PackedFragmentBuffer.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace build
{

/**
 * \brief This class implements the sorting for single-ended data:
 **/
class NotAFilter
{
public:
    template<typename InputIteratorT, typename InsertIteratorT>
    void filterInput(
        const PackedFragmentBuffer &fragments,
        InputIteratorT duplicatesBegin,
        InputIteratorT duplicatesEnd,
        BuildStats &buildStats,
        const unsigned binIndex,
        InsertIteratorT results)
    {
        BOOST_FOREACH(const typename InputIteratorT::value_type &idx, std::make_pair(duplicatesBegin, duplicatesEnd))
        {
            const io::FragmentAccessor &fragment = fragments.getFragment(idx);
            buildStats.incrementUniqueFragments(binIndex, fragment.barcode_);
            buildStats.incrementTotalFragments(binIndex, fragment.barcode_);
            results++ = PackedFragmentBuffer::Index(idx, fragment);
        }
    }
};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_NOT_A_FILTER_HH
