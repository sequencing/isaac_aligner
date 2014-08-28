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
 ** \file ReferenceSorter.hh
 **
 ** Top level component to produce a sorted reference.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_REFERENCE_SORTER_HH
#define iSAAC_REFERENCE_REFERENCE_SORTER_HH

#include <string>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>

#include "oligo/Kmer.hh"
#include "reference/ReferenceKmer.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace reference
{

template <typename KmerT>
class ReferenceSorter: boost::noncopyable
{
public:
    ReferenceSorter(
        const unsigned int maskWidth,
        const unsigned mask,
        const boost::filesystem::path &genomeFile,
        const boost::filesystem::path &genomeNeighborsFile,
        const boost::filesystem::path &outputFile,
        const unsigned repeatThreshold);
    void run();
private:
    const unsigned repeatThreshold_;
    const unsigned int maskWidth_;
    const unsigned mask_;

    // the mask highlight bits in original kmer (ABCD)
    const KmerT msbMask_;
    // the mask value in the original kmer (ABCD) shifted to the topmost position
    const KmerT maskBits_;

    // the mask highlight bits as if the permutated kmer having them all set was unpermutated back into ABCD
    const KmerT unpermutatedMsbMask_;
    // the mask as if the permutated kmer containing it was unpermutated back into ABCD
    const KmerT unpermutatedMaskBits_;

    const boost::filesystem::path genomeFile_;
    const boost::filesystem::path genomeNeighborsFile_;

    const boost::filesystem::path outputFile_;
    std::vector<ReferenceKmer<KmerT> > reference_;


    unsigned long loadReference(
        std::vector<unsigned long> &contigOffsets);
    void sortReference();
    void saveReference(
        const std::vector<unsigned long> &contigOffsets,
        const std::vector<bool> &neighbors,
        const unsigned long genomeLength);
    void addToReference(const KmerT kmer, const ReferencePosition &referencePosition);
};

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_REFERNECE_SORTER_HH
