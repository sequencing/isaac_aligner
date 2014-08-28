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
 ** \file NeighborsFinder.hh
 **
 ** \brief Top level component to find neighbors.
 **
 ** \author Come Raczy
 **/

#ifndef ISAAC_REFERENCE_NEIGHBORS_FINDER_HH
#define ISAAC_REFERENCE_NEIGHBORS_FINDER_HH

#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>

#include "oligo/Kmer.hh"
#include "oligo/Permutate.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace reference
{

template <typename KmerT>
class NeighborsFinder: boost::noncopyable
{
public:
    struct AnnotatedKmer
    {
        AnnotatedKmer(const KmerT kmer, bool hasNghbrs) : value(kmer), hasNeighbors(hasNghbrs) {}
        KmerT value;
        bool hasNeighbors;
        bool operator<(const AnnotatedKmer &rhs) const {return this->value < rhs.value;}
        void setHasNeighbors(){hasNeighbors = true;}
    }__attribute__ ((packed)); // this significantly reduces memory requirement for finder especially considering the fact
                               // that parallel sort needs twice the memory for processing.

    typedef std::vector<AnnotatedKmer> KmerList;
    NeighborsFinder(
        const bool parallelSort,
        const boost::filesystem::path &inputFile,
        const boost::filesystem::path &outputDirectory,
        const boost::filesystem::path &outputFile,
        const boost::filesystem::path &tempFile,
        const unsigned jobs);
    void run() const;
    static void findNeighbors(KmerList &kmerList, unsigned jobs);
    /**
     ** \brief Count the non-equal neighbors within Hamming distance of neighborhoodWidth
     **
     ** The whole block must share the same prefix as the given kmer
     **/
    static void markNeighbors(
        const typename KmerList::iterator blockBegin,
        const typename KmerList::const_iterator blockEnd);
private:
    const bool parallelSort_;
    const boost::filesystem::path inputFile_;
    const boost::filesystem::path outputDirectory_;
    const boost::filesystem::path outputFile_;
    const boost::filesystem::path tempFile_;
    const unsigned jobs_;
    static const unsigned neighborhoodWidth = 4;

    void generateNeighbors(const SortedReferenceMetadata &sortedReferenceMetadata) const;
    void storeNeighborKmers(const KmerList &kmerList) const;
    void updateSortedReference(SortedReferenceMetadata::MaskFiles &maskFileList) const;
    static void findNeighborsParallel(const typename KmerList::iterator kmerListBegin, const typename KmerList::iterator kmerListEnd);
    KmerList getKmerList(const SortedReferenceMetadata &sortedReferenceMetadata) const;
};

} // namespace reference
} //namespace isaac

#endif // #ifndef ISAAC_REFERENCE_NEIGHBORS_FINDER_HH
