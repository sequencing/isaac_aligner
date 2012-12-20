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

/**
 * \brief defines match for the purpose of the alignment.
 */
inline bool isMatch(const char readBase, const char referenceBase)
{
    return readBase == referenceBase || referenceBase == 'N' || readBase == 'n';
}

/**
 * \brief moves sequenceBegin to the first position followed by CONSECUTIVE_MATCHES_MAX matches
 *
 * \return distance moved
 */
template <unsigned CONSECUTIVE_MATCHES_MIN, typename SequenceIteratorT, typename ReferenceIteratorT, typename BaseExtractor>
unsigned clipMismatches(
    SequenceIteratorT sequenceBegin, const SequenceIteratorT sequenceEnd,
    ReferenceIteratorT referenceBegin, ReferenceIteratorT referenceEnd,
    BaseExtractor baseExtractor)
{
    unsigned matchesInARow = 0;
    unsigned ret = 0;
    while (sequenceEnd != sequenceBegin && referenceBegin != referenceEnd && CONSECUTIVE_MATCHES_MIN > matchesInARow)
    {
        if (isMatch(baseExtractor(*sequenceBegin), *referenceBegin))
        {
            ++matchesInARow;
        }
        else
        {
            matchesInARow = 0;
        }
        ++sequenceBegin;
        ++referenceBegin;
        ++ret;
    }

    return (CONSECUTIVE_MATCHES_MIN == matchesInARow) ? ret - matchesInARow: 0;
}



} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_ALIGNMENT_HH
