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
 ** \file SingleEndFilter.hh
 **
 ** Template for general approach to filtering ends of duplicate pairs.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_SINGLE_END_FILTER_HH
#define iSAAC_BUILD_SINGLE_END_FILTER_HH

#include "build/BuildStats.hh"
#include "build/PackedFragmentBuffer.hh"
#include "common/Debug.hh"
#include "io/FragmentIndex.hh"

namespace isaac
{
namespace build
{

/**
 * \brief This class implements the sorting for single-ended data:
 **/
class SingleEndFilter
{
public:
    template<typename InsertIteratorT>
    void filterInput(
        const PackedFragmentBuffer &fragments,
        std::vector<io::SeFragmentIndex>::iterator duplicatesBegin,
        std::vector<io::SeFragmentIndex>::iterator duplicatesEnd,
        BuildStats &buildStats,
        const unsigned binIndex,
        InsertIteratorT results)
    {
        BOOST_FOREACH(const io::SeFragmentIndex &idx, std::make_pair(duplicatesBegin, duplicatesEnd))
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

#endif // #ifndef iSAAC_BUILD_SINGLE_END_FILTER_HH
