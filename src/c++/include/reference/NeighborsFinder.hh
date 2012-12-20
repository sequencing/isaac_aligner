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
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace reference
{

typedef SortedReferenceXml::MaskFile MaskFile;

class NeighborsFinder: boost::noncopyable
{
public:
    struct AnnotatedKmer
    {
        AnnotatedKmer(const oligo::Kmer kmer, bool hasNghbrs) : value(kmer), hasNeighbors(hasNghbrs) {}
        oligo::Kmer value;
        bool hasNeighbors;
        bool operator<(const AnnotatedKmer &rhs) const {return this->value < rhs.value;}
        void setHasNeighbors(){hasNeighbors = true;}
    }__attribute__ ((packed)); // this significantly reduces memory requirement for finder especially considering the fact
                               // that standard parallel sort needs twice the memory for processing.
    typedef std::vector<AnnotatedKmer> KmerList;
    NeighborsFinder(
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
        const KmerList::iterator blockBegin,
        const KmerList::const_iterator blockEnd);
private:
    const boost::filesystem::path inputFile_;
    const boost::filesystem::path outputDirectory_;
    const boost::filesystem::path outputFile_;
    const boost::filesystem::path tempFile_;
    const unsigned jobs_;
    static const unsigned neighborhoodWidth = 4;

    void generateNeighbors(const SortedReferenceXml &sortedReferenceXml) const;
    void storeNeighborKmers(const KmerList &kmerList) const;
    void updateSortedReference(std::vector<MaskFile> &maskFileList) const;
    static void findNeighborsParallel(const KmerList::iterator kmerListBegin, const KmerList::iterator kmerListEnd);
};

} // namespace reference
} //namespace isaac

#endif // #ifndef ISAAC_REFERENCE_NEIGHBORS_FINDER_HH
