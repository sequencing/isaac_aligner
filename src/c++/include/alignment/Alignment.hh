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
 ** \file Alignment.hh
 **
 ** \brief Basic alignment constants and utilities.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_ALIGNMENT_HH
#define iSAAC_ALIGNMENT_ALIGNMENT_HH

#include <string>
#include <vector>
#include <stdint.h>

#include "alignment/Read.hh"
#include "alignment/Quality.hh"

namespace isaac
{
namespace alignment
{

// Short seeds require extra evidence such as absence of seed neighbors in the reference or having more than on non-overlapping seed
static const unsigned WEAK_SEED_LENGTH = 32;
// Long seeds are good enough evidence for non-ambiguous anchoring
static const unsigned STRONG_SEED_LENGTH = 64;


/**
 * \brief defines match for the purpose of the alignment.
 */
inline bool isMatch(const char readBase, const char referenceBase)
{
    return readBase == 'n' || (readBase == referenceBase && referenceBase != 'N');
}

/**
 * \brief moves sequenceBegin to the first position followed by CONSECUTIVE_MATCHES_MAX matches
 *
 * \return pair(distance moved, edit distance adjustment). edit distance adjustment equals all mismatches that have been clipped away.
 *         Note that N is considered to be an edit distance mismatch in this case.
 */
template <unsigned CONSECUTIVE_MATCHES_MIN, typename SequenceIteratorT, typename ReferenceIteratorT, typename BaseExtractor>
std::pair<unsigned, unsigned> clipMismatches(
    SequenceIteratorT sequenceBegin, const SequenceIteratorT sequenceEnd,
    ReferenceIteratorT referenceBegin, ReferenceIteratorT referenceEnd,
    BaseExtractor baseExtractor)
{
    unsigned matchesInARow = 0;
    unsigned ediDistanceMismatches = 0;
    // The number of mismatches that are not part of the sequence that gets clipped
    unsigned ediDistanceMismatchesUnclipped = 0;
    unsigned ret = 0;
    while (sequenceEnd != sequenceBegin && referenceBegin != referenceEnd && CONSECUTIVE_MATCHES_MIN > matchesInARow)
    {
        char sequenceBase = baseExtractor(*sequenceBegin);
        if (isMatch(sequenceBase, *referenceBegin))
        {
            ++matchesInARow;
            ediDistanceMismatchesUnclipped += (sequenceBase != *referenceBegin);
        }
        else
        {
            matchesInARow = 0;
            ediDistanceMismatchesUnclipped = 0;
        }
        ediDistanceMismatches += (sequenceBase != *referenceBegin);
        ++sequenceBegin;
        ++referenceBegin;
        ++ret;
    }

    return (CONSECUTIVE_MATCHES_MIN == matchesInARow) ?
        std::make_pair(ret - matchesInARow, ediDistanceMismatches - ediDistanceMismatchesUnclipped):
        std::make_pair(0U,0U);
}

template <typename SequenceIteratorT, typename BaseExtractor>
unsigned countMatches(
    SequenceIteratorT sequenceBegin,
    const SequenceIteratorT sequenceEnd,
    std::vector<char>::const_iterator referenceBegin,
    const std::vector<char>::const_iterator referenceEnd,
    BaseExtractor baseExtractor)
{
    unsigned ret = 0;
    for (;sequenceEnd != sequenceBegin && referenceEnd != referenceBegin;
        ++sequenceBegin, ++referenceBegin)
    {
        ret += isMatch(baseExtractor(*sequenceBegin), *referenceBegin);
    }
    return ret;
}

template <typename SequenceIteratorT>
unsigned countMatches(
    const SequenceIteratorT sequenceBegin,
    const SequenceIteratorT sequenceEnd,
    const std::vector<char>::const_iterator referenceBegin,
    const std::vector<char>::const_iterator referenceEnd)
{
    return countMatches(sequenceBegin, sequenceEnd,
                        referenceBegin, referenceEnd,
                        &boost::cref<typename std::iterator_traits<SequenceIteratorT>::value_type>);
}

template <typename SequenceIteratorT, typename BaseExtractor>
unsigned countMismatches(
    SequenceIteratorT sequenceBegin,
    const SequenceIteratorT sequenceEnd,
    std::vector<char>::const_iterator referenceBegin,
    const std::vector<char>::const_iterator referenceEnd,
    BaseExtractor baseExtractor)
{
    unsigned ret = 0;
    for (;sequenceEnd != sequenceBegin && referenceEnd != referenceBegin;
        ++sequenceBegin, ++referenceBegin)
    {
        ret += !isMatch(baseExtractor(*sequenceBegin), *referenceBegin);
    }
    return ret;
}

template <typename SequenceIteratorT>
unsigned countMismatches(
    const SequenceIteratorT sequenceBegin,
    const SequenceIteratorT sequenceEnd,
    const std::vector<char>::const_iterator referenceBegin,
    const std::vector<char>::const_iterator referenceEnd)
{
    return countMismatches(sequenceBegin, sequenceEnd,
                           referenceBegin, referenceEnd,
                           &boost::cref<typename std::iterator_traits<SequenceIteratorT>::value_type>);
}


template <typename SequenceIteratorT, typename BaseExtractor>
unsigned countMismatches(
    const SequenceIteratorT basesIterator,
    const std::vector<char>::const_iterator referenceBegin,
    const std::vector<char>::const_iterator referenceEnd,
    unsigned length,
    BaseExtractor baseExtractor)
{
    return countMismatches(basesIterator, basesIterator + length,
                           referenceBegin, referenceEnd, baseExtractor);
}

} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_ALIGNMENT_HH
