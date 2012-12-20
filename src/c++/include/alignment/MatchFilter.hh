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
 ** \file MatchFilter.hh
 **
 ** \brief Filtering mechanism for the MatchFinder, based on the number of
 ** mismatches in each block.
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_FILTER_HH
#define iSAAC_ALIGNMENT_MATCH_FILTER_HH

#include <vector>
#include <string>
#include <boost/noncopyable.hpp>

namespace isaac
{
namespace alignment
{

/**
 ** \brief Filtering mechanism for matches, based on the number of mismatches in
 ** each block
 **
 ** The aligner allows up to two mismatches. To do so, it splits the k-mer into
 ** 4 blocks, ABCD, finds an exact match on the first two blocks (AB) and counts
 ** the mismatches in the other two blocks (CD). The valid number of mismatches
 ** in C and D respectively are (0, 0), (0, 1), (0, 2), (1, 0), (1, 1) and (2,
 ** 0). To dind all possible matches with up-to two mismatches, the aligner
 ** repeats the procedure over six permuatrions of the kmer: ABCD, BCDA, CDAB,
 ** ACBD, BDAC and ADBC (that list is defined in "oligo/Permutations.hh".
 **
 ** The purpose of this filter is to prevent the aligner to secord repeatedly
 ** the matches that were already found in the previous permutations. for
 ** instance, the exact matches need to be produced only when matching the
 ** original k-mer ABCD. Similarly, when processing BCDA, all matches with an
 ** exact match on A and up to two mismatches on D have already bee found when
 ** processing ABCD and should be filtered out.
 **/
class MatchFilter: boost::noncopyable
{
public:
    MatchFilter(const std::string permutation) : use_(getUse(permutation)) {}
    bool use(unsigned int mismatchCount1, unsigned int mismatchCount2)
    {
        return
            (mismatchCount1 + mismatchCount2) <= 2 &&
            use_[(mismatchCount1 << 2) | mismatchCount2];
    }
private:
    MatchFilter();
    const std::vector<bool> use_;
    std::vector<bool> getUse(const std::string permutation);
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_FILTER_HH
