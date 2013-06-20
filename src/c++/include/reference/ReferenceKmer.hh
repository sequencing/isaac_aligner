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

template <typename KmerT>
class ReferenceKmer
{
public:
    KmerT first;
    ReferencePosition::value_type second;
    /// Encodes a contig Id and a position into a ReferencePosition
    ReferenceKmer(
        const KmerT &kmer = 0,
        const ReferencePosition &referencePosition = ReferencePosition(0))
        : first(kmer), second(referencePosition.getValue())
    {
    }
    KmerT getKmer() const {return first;}
    const ReferencePosition getTranslatedPosition(
        const std::vector<unsigned> &contigTranslationTable) const {return ReferencePosition(second).translateContig(contigTranslationTable);}
    const ReferencePosition getReferencePosition() const {return ReferencePosition(second);}
    void setKmer(const KmerT kmer) {first = kmer;}
    void setNeighbors(bool set){second = ReferencePosition(second).setNeighbors(set).getValue();}
    void setNeighbors(){second = ReferencePosition(second).setNeighbors(true).getValue();}
    bool hasNoNeighbors() const {return !ReferencePosition(second).hasNeighbors();}
} __attribute__ ((packed));

BOOST_STATIC_ASSERT(sizeof(ReferenceKmer<oligo::KmerType>) == (sizeof(oligo::KmerType) + sizeof(ReferencePosition)));
BOOST_STATIC_ASSERT(sizeof(ReferenceKmer<oligo::LongKmerType>) == (sizeof(oligo::LongKmerType) + sizeof(ReferencePosition)));

template <typename KmerT>
inline std::ostream &operator<<(std::ostream &os, const ReferenceKmer<KmerT> &rk)
{
    return os << "ReferenceKmer(" << isaac::oligo::bases(rk.getKmer()) << "," << rk.getReferencePosition() << ")";
}


template <typename KmerT>
inline bool compareKmer(const ReferenceKmer<KmerT> &lhs, const ReferenceKmer<KmerT> &rhs)
{
    return lhs.getKmer() < rhs.getKmer();
}

template <typename KmerT>
inline bool comparePosition(const ReferenceKmer<KmerT> &lhs, const ReferenceKmer<KmerT> &rhs)
{
    return lhs.getReferencePosition() < rhs.getReferencePosition();
}

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_REFERENCE_KMER_HH
