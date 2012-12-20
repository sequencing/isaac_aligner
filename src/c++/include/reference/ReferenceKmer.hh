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
 ** \file ReferenceKmer.hh
 **
 ** Representation of a k-mer at a given position in a reference genome.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_REFERENCE_KMER_HH
#define iSAAC_REFERENCE_REFERENCE_KMER_HH

#include <utility>

#include "oligo/Kmer.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace reference
{

class ReferenceKmer: public std::pair<oligo::Kmer, ReferencePosition>
{
public:
    /// Encodes a contig Id and a position into a ReferencePosition
    ReferenceKmer(
        const oligo::Kmer &kmer = 0,
        const ReferencePosition &referencePosition = ReferencePosition(0))
        : std::pair<oligo::Kmer, ReferencePosition>(kmer, referencePosition)
    {
    }
    oligo::Kmer getKmer() const {return first;}
    const ReferencePosition getTranslatedPosition(
        const std::vector<unsigned> &contigTranslationTable) const {return second.translateContig(contigTranslationTable);}
    const ReferencePosition &getReferencePosition() const {return second;}
    ReferencePosition &getReferencePosition() {return second;}
    void setKmer(oligo::Kmer kmer) {first = kmer;}
    void setNeighbors(){second.setNeighbors(true);}
};

BOOST_STATIC_ASSERT(16 == sizeof(ReferenceKmer));

inline std::ostream &operator<<(std::ostream &os, const ReferenceKmer &rk)
{
    return os << "ReferenceKmer(" << rk.getKmer() << "," << rk.getReferencePosition() << ")";
}


inline bool compareKmer(const ReferenceKmer &lhs, const ReferenceKmer &rhs)
{
    return lhs.getKmer() < rhs.getKmer();
}
inline bool comparePosition(const ReferenceKmer &lhs, const ReferenceKmer &rhs)
{
    return lhs.getReferencePosition() < rhs.getReferencePosition();
}

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_REFERENCE_KMER_HH
