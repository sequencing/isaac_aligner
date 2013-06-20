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
 ** \file BandedSmithWaterman.hh
 **
 ** \brief SSE2 implementation of a banded smith waterman
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_BANDED_SMITH_WATERMAN_HH
#define iSAAC_ALIGNMENT_BANDED_SMITH_WATERMAN_HH

#include <string>

#include <boost/noncopyable.hpp>

#include "alignment/Cigar.hh"

namespace isaac
{
namespace alignment
{
    
/** 
 ** \brief global optimization for alignments with a maximum gap size
 ** 
 ** Assumes that the scores fit into a short integer (2 bytes).
 ** Requires SSE4 for "v8hi __builtin_ia32_pmaxsw128 (v8hi, v8hi)".
 ** _mm_max_pi16 is probably an mmx instruction
 ** 
 ** This implementation is limited to a maximum gap width of 15 bases
 ** (+/-7 bands) so that all operations can be done over two registers.
 **
 ** The registers are aligned to the database.
 **
 ** Note: this is non-copyable because of the dynamically-allocated internal
 ** buffer.
 ** 
 **/
class BandedSmithWaterman: boost::noncopyable
{
public:
    /**
     * \brief Initialize the optimizer with specific scores and width
     *
     * \param matchScore - Expected to be positive. The higher the value, the more likely matches are chosen
     * \param mismatchScore - Expected to be negative. The lower the value, the less likely the mismatches are chosen
     * \param gapOpenScore - Expected to be positive. The higher the value, the less likely the gaps are opened
     * \param gapOpenScore - Expected to be positive. The higher the value, the less likely the gaps are extended
     */
    BandedSmithWaterman(
        int matchScore, int mismatchScore, int gapOpenScore,
        int gapExtendScore, int maxReadLength);
    /// \brief delete the pre-allocated re-usable buffer
    ~BandedSmithWaterman();
    /**
     ** \brief align the query to the database and store the descriptor of the best match.
     **
     ** This is the translation of implementation from ELAND into SIMD, but
     ** extending the widest gap size to 7 (instead of 5) because this comes for
     ** free with the size of the SIMD registers.
     **
     ** Note: this operation is not 'const' because it uses a pre-allocated internal buffer.
     **/
    unsigned align(
        const std::vector<char> &query,
        const std::vector<char>::const_iterator databaseBegin,
        const std::vector<char>::const_iterator databaseEnd,
        Cigar &cigar) const;

    unsigned align(
        const std::vector<char>::const_iterator queryBegin,
        const std::vector<char>::const_iterator queryEnd,
        const std::vector<char>::const_iterator databaseBegin,
        const std::vector<char>::const_iterator databaseEnd,
        Cigar &cigar) const;

    // the widest gap-size handled by this implementation
    static const unsigned WIDEST_GAP_SIZE = 16;
    // if we know there are no reference matching kmers within cutoffDistance,
    // there is no point to do the gapped alignment.
    static const unsigned distanceCutoff = 7;
    // Assumedly no reason to do gapped alignment if total mismatch count is 5 or less
    static const unsigned mismatchesCutoff = 5;
private:
    const int matchScore_;
    const int mismatchScore_; 
    const int gapOpenScore_;
    const int gapExtendScore_;
    const int maxReadLength_;
    const short initialValue_; // minimal usable value to initialize the matrices
    typedef unsigned short ScoreType;
    static const unsigned int registerLength_ = 16 / sizeof(ScoreType);
    char *T_;
};  

} // namespace alignment
} // namespace isaac
    
#endif // #ifndef iSAAC_ALIGNMENT_BANDED_SMITH_WATERMAN_HH
