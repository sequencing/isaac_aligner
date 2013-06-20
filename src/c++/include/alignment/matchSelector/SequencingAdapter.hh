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
 ** \file SequencingAdapter.hh
 **
 ** \brief Helper class for verifying adapter sequence matches
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_SEQUENCING_ADAPTER_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_SEQUENCING_ADAPTER_HH

#include <string>
#include <vector>

#include "flowcell/SequencingAdapterMetadata.hh"
#include "oligo/KmerGenerator.hpp"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

/**
 * \brief Helper class for verifying adapter sequence matches
 */
class SequencingAdapter
{
    static const unsigned adapterMatchBasesMin_ = 5;
    static const char UNINITIALIZED_POSITION = -1;
    static const char NON_UNIQUE_KMER_POSITION = -2;

    flowcell::SequencingAdapterMetadata adapterMetadata_;

    std::vector<char> kmerPositions_;
public:
    SequencingAdapter(const flowcell::SequencingAdapterMetadata &adapterMetadata);

    const std::pair<std::vector<char>::const_iterator, std::vector<char>::const_iterator>
    getMatchRange(
        const std::vector<char>::const_iterator sequenceBegin,
        const std::vector<char>::const_iterator sequenceEnd,
        const std::vector<char>::const_iterator mismatchBase) const;

    /**
     * \brief Unbounded adapters can be only found on the strand which they match.
     *        fixed-length adapters can be found on any strand in the order in which
     *        they appear in the list of adapters for the sample prep
     */
    bool isStrandCompatible(const bool reverse) const
    {
        return !adapterMetadata_.isUnbounded() || reverse == adapterMetadata_.isReverse();
    }

private:
    static bool isGoodPosition(char pos) {return 0 <= pos;}
};

typedef std::vector<SequencingAdapter> SequencingAdapterList;

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_SEQUENCING_ADAPTER_HH
